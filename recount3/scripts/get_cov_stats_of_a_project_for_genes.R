#!/usr/bin/env Rscript
#
# QC項目の説明：https://rna.recount.bio/docs/quality-check-fields.html
# 「1. プロジェクト情報の取得」だけを実行すればデータをダウンロードしてキャッシュに蓄積できる、これを並列でやるとエラー出るから直列でする

suppressPackageStartupMessages({
  library(argparse)
  library(recount3)
  library(megadepth)
  library(rtracklayer)
  library(GenomicRanges)
  library(dplyr)
})


parser <- ArgumentParser()
parser$add_argument("--gtf", type="character", required=TRUE, help="Path to GTF annotation file")
parser$add_argument("--gene_ids", type="character", default=NULL, help="Path to file with gene_ids to process (one per line). If omitted, all genes in GTF are used.")
parser$add_argument("--project", type="character", required=TRUE, help="human_projectsの一列目で絞るためのプロジェクト名")
parser$add_argument("--seed", type="integer", default=123,
                    help="Random seed for sampling [default %(default)s]")
parser$add_argument("--n_sample", type="integer", default=10,
                    help="Number of files to sample [default %(default)s]")
parser$add_argument("--threshold", type="integer", default=10,
                    help="Row count threshold to trigger sampling [default %(default)s]")
parser$add_argument("--output", type="character", required=TRUE)
args <- parser$parse_args()


# 1. GTFからエクソン領域を抽出し、遺伝子ごとにreduceして重複をなくす
message("1: Loading GTF and preparing exon regions")
gtf <- import(args$gtf)
exons <- gtf[gtf$type == "exon"]

# gene_idリストが指定されていればフィルタリング
if (!is.null(args$gene_ids)) {
  target_ids <- readLines(args$gene_ids)
  target_ids <- target_ids[nchar(trimws(target_ids)) > 0]  # 空行除去
  exons <- exons[exons$gene_id %in% target_ids]
  message(paste("  Filtering to", length(target_ids), "gene_ids ->", length(unique(exons$gene_id)), "found in GTF"))
}

if (length(exons) == 0) {
  stop("No exons found after filtering. Check --gtf and --gene_ids arguments.")
}

# 遺伝子ごとにエクソンをreduceして重複領域をマージ（二重カウント防止）
exon_by_gene <- split(exons, exons$gene_id)
reduced_by_gene <- GenomicRanges::reduce(exon_by_gene)
flat_exons <- unlist(reduced_by_gene)

# 座標→gene_idの対応テーブル（行順ではなく座標でjoinするため）
# get_coverage()はBED入力の行順と異なる座標順でGRangesを返すことがある
coord_to_gene <- data.frame(
  seqnames = as.character(seqnames(flat_exons)),
  start    = start(flat_exons),  # GRanges 1-based
  end      = end(flat_exons),
  gene_id  = names(flat_exons),
  stringsAsFactors = FALSE
)

# megadepth用の一時BEDファイルを作成（0-based half-open）
# 0-basedにした方が正確だが後でjoinするときに戻さないといけないのでそのままで
tmp_bed <- tempfile(fileext = ".bed")
exon_bed_df <- data.frame(
  chrom  = coord_to_gene$seqnames,
  start  = coord_to_gene$start,
  end    = coord_to_gene$end,
  name   = coord_to_gene$gene_id,
  score  = 0,
  strand = "."
)
write.table(exon_bed_df, tmp_bed, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
message(paste("  Exon regions:", nrow(exon_bed_df), "/ Genes:", length(unique(coord_to_gene$gene_id))))


# 2. プロジェクト情報の取得
message("2: Fetching project metadata")
human_projects <- available_projects()
dataset_info <- filter(human_projects, project == args$project)
rse <- create_rse(dataset_info)
metadata_df <- as.data.frame(colData(rse)) %>%
  select(BigWigURL, recount_qc.bc_auc.unique_reads_all_bases)


# 3. サンプリング
message("3: Sampling")
all_indices <- seq_len(nrow(metadata_df))

if (length(all_indices) > args$threshold) {
  set.seed(args$seed)
  sample_n <- min(length(all_indices), args$n_sample)
  sampled_indices <- sample(all_indices, size = sample_n)
  target_data <- metadata_df[sampled_indices, ]
  message(paste("  Sampling", sample_n, "files from", length(all_indices)))
} else {
  target_data <- metadata_df
  message(paste("  Using all", length(all_indices), "files (below threshold)"))
}


# 4. bigWigからエクソンカバレッジを取得し、遺伝子ごとに集計
message("4: Computing coverage per gene")
results_list <- list()

for (i in seq_len(nrow(target_data))) {
  bw_url      <- target_data$BigWigURL[i]
  norm_factor <- target_data$recount_qc.bc_auc.unique_reads_all_bases[i]

  message(paste("  Processing:", bw_url))

  # エクソンごとのカバレッジ（座標でcoord_to_geneとjoinしてgene_idを付ける）
  res    <- get_coverage(bw_url, op = "sum", annotation = tmp_bed)
  res_df <- as.data.frame(res) %>%
    mutate(seqnames = as.character(seqnames)) %>%
    left_join(coord_to_gene, by = c("seqnames", "start", "end"))

  # 遺伝子ごとにエクソンのscoreを合算してからAUC正規化
  gene_df <- res_df %>%
    group_by(gene_id) %>%
    summarise(score = sum(score), .groups = "drop") %>%
    mutate(
      score_normalized = (score / norm_factor) * 1e6,
      bw_url = bw_url
    )

  results_list[[bw_url]] <- gene_df
}


# 5. 全データの結合と保存
message("5: Saving results")
final_df <- bind_rows(results_list)

out_dir <- dirname(args$output)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

write.table(final_df, args$output, sep = "\t", row.names = FALSE, quote = FALSE)
message(paste("Success! Saved to:", args$output))
