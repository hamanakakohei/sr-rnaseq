#!/usr/bin/env Rscript
#
# QC項目の説明：https://rna.recount.bio/docs/quality-check-fields.html
# 「1. プロジェクト情報の取得」だけを実行すればデータをダウンロードしてキャッシュに蓄積できる、これを並列でやるとエラー出るから直列でする

suppressPackageStartupMessages({
  library(argparse)
  library(recount3)
  library(megadepth)
  library(dplyr)
})


parser <- ArgumentParser()
parser$add_argument("--op", type="character", default="sum", help="Operation for get_coverage: sum, mean, max, min [default %(default)s]")
parser$add_argument("--annotation", type="character", required=TRUE, help="Path to the annotation file (BED)")
parser$add_argument("--project", type="character", required=TRUE, help="human_projectsの一列目で絞るためのプロジェクト名")
parser$add_argument("--seed", type="integer", default=123,
                    help="Random seed for sampling [default %(default)s]")
parser$add_argument("--n_sample", type="integer", default=10,
                    help="Number of files to sample [default %(default)s]")
parser$add_argument("--threshold", type="integer", default=10,
                    help="Row count threshold to trigger sampling [default %(default)s]")
parser$add_argument("--output", type="character", required=TRUE)
args <- parser$parse_args()


# 1. プロジェクト情報の取得
human_projects <- available_projects()
dataset_info <- filter(human_projects, project == args$project)
rse <- create_rse(dataset_info)
metadata_df <- as.data.frame(colData(rse)) %>%
  select(BigWigURL, recount_qc.bc_auc.unique_reads_all_bases)


# 2. サンプリング
# 行数（ファイル数）が閾値を超えている場合のみランダムサンプリング
message("2")
all_indices <- seq_len(nrow(metadata_df))

if (length(all_indices) > args$threshold) {
  set.seed(args$seed)
  sample_n <- min(length(all_indices), args$n_sample)
  sampled_indices <- sample(all_indices, size = sample_n)
  target_data <- metadata_df[sampled_indices, ]
  message(paste("Sampling", sample_n, "files from", length(all_indices)))
} else {
  target_data <- metadata_df
  message(paste("Using all", length(all_indices), "files (below threshold)"))
}


# 3. データの取得と結合
# 結果をリストに格納していく
message("3")
results_list <- list()

for(i in 1:nrow(target_data)){
  bw_url <- target_data$BigWigURL[i]
  norm_factor <- target_data$recount_qc.bc_auc.unique_reads_all_bases[i]

  message(paste("Processing:", bw_url))

  res <- get_coverage(bw_url, op = args$op, annotation = args$annotation)

  res_df <- as.data.frame(res) %>%
    mutate(
      score = (score / norm_factor) * 1e6,
      bw_url = bw_url
    )

  results_list[[bw_url]] <- res_df
}


# 4. 全データの結合と保存
# bind_rowsでリスト内の全データフレームを上下に結合
# ディレクトリが存在しない場合は作成
message("4")
final_df <- bind_rows(results_list)

out_dir <- dirname(args$output)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

write.table(final_df, args$output, sep = "\t", row.names = FALSE, quote = FALSE)
message(paste("Success! Saved to:", args$output))
