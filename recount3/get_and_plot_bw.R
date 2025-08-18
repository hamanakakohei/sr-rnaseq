#!/usr/bin/env Rscript

# 参考：
# https://bioconductor.org/packages/release/bioc/vignettes/recount3/inst/doc/recount3-quickstart.html
# https://bioconductor.org/packages/3.21/bioc/vignettes/derfinder/inst/doc/derfinder-users-guide.html
# https://bioconductor.org/packages/3.21/bioc/vignettes/derfinderPlot/inst/doc/derfinderPlot.html

library("recount3")
suppressPackageStartupMessages(library("derfinder"))
#library("derfinderData")
library("derfinderPlot")
suppressPackageStartupMessages(library("TxDb.Hsapiens.UCSC.hg38.knownGene"))
library(bumphunter)
library(GenomicRanges)


# --- 引数 ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript get_and_plot_bw.R <SRP_ID> <bw_and_plot_dir> <seqnames> <starts> <ends> [bw_overwrite=FALSE]")
  # 例えば：
  # Rscript ./get_and_plot_bw.R \
  #   SRP10000 \
  #   results \
  #   "chr1,chr1" \
  #   "1000,2000" \
  #   "3000,4000" \
  #   FALSE
}

SRP_ID   <- args[1]
outdir   <- args[2]
seqnames <- strsplit(args[3], ",")[[1]]
starts   <- as.integer(strsplit(args[4], ",")[[1]])
ends     <- as.integer(strsplit(args[5], ",")[[1]])
overwrite <- ifelse(length(args) >= 6 && args[6] == "TRUE", TRUE, FALSE)

print(args[3])
print(seqnames)
print(args[4])
print(starts)
print(args[5])
print(ends)

dir.create(outdir, showWarnings = FALSE)


# --- プロジェクトを指定して、rseとbw URLベクターを作る ---
human_projects <- available_projects()
project_info <- subset(
    human_projects,
    project == SRP_ID & project_type == "data_sources"
) 
rse_gene <- create_rse(project_info)

files <- setNames(rse_gene$BigWigURL, rse_gene$external_id)


# --- 通信の問題で、URLに直接fullCoverage()を使えないとき---
# 保存先ディレクトリ
outdir <- "bw"
dir.create(outdir, showWarnings = FALSE)

# ダウンロード
mapply(function(url, sample) {
  outfile <- file.path(outdir, basename(url))
  if (!file.exists(outfile) || overwrite) {
    #cmd <- sprintf("aria2c -x16 -s16 -o %s -d %s %s", basename(url), outdir, url)
    cmd <- sprintf("wget -O %s %s", outfile, url)
    #cmd <- sprintf("curl -L -o %s %s", outfile, url)
    system(cmd)
  }
  outfile
}, files, names(files)) -> local_paths

files <- setNames(local_paths, names(files))


# --- もろもろの準備 ---
# regionsを作る
regions <- GRanges(
  seqnames = seqnames,
  ranges   = IRanges(start = starts, end = ends),
  strand   = "*"
)

# regionCoverage を作る
fullCov <- fullCoverage(
    files = files, 
    chrs = unique(seqnames),
    totalMapped = unname(sapply(files, getTotalMapped))
    #targetSize = 5e7
)

regionCov <- getRegionCoverage(fullCov, regions)

# nearestAnnotation を作る
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes <- annotateTranscripts(txdb = TxDb.Hsapiens.UCSC.hg38.knownGene)
nearestAnnotation <- matchGenes(x = regions, subject = genes)

# annotatedRegions を作る
annoRegs <- annotateRegions(regions, genomicState$fullGenome)


# --- groupInfoを作る ---
print(rse_gene$sra.experiment_title)
print(rse_gene$sra.sample_title)

## 例１：sra.sample_attributesから自動で作る
#rse_gene_sraInfo <- expand_sra_attributes(rse_gene)
#rse_gene_sraInfo = colData(rse_gene_sraInfo)[, (ncol(colData(rse_gene))+1):ncol(colData(rse_gene_sraInfo))]
#group <- apply(as.data.frame(rse_gene_sraInfo), 1, paste, collapse = ";")
#group <- factor(group)

## 例2：sra.experiment_title や sra.sample_title を見てマニュアル編集
group <- factor(c(
  "WT", "WT", "WT",
  "KD", "KD", "KD"
))


# --- 各 region ごとに画像保存 ---
for (i in seq_along(regions)) {
  png(file.path(outdir, paste0("region_", i, ".png")), width=1200, height=800)
  plotRegionCoverage(
    regions = regions,
    regionCoverage = regionCov,
    groupInfo = group,
    nearestAnnotation = nearestAnnotation,
    annotatedRegions = annoRegs,
    ask = FALSE, verbose = FALSE,
    txdb = txdb, scalefac = 1,
    whichRegions = i
  )
  dev.off()
}
