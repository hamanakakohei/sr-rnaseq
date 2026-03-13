

## Load recount3 R package
library("recount3")

## Find all available human projects
human_projects <- available_projects()

## Find the project you are interested in,
## here we use SRP009615 as an example
proj_info <- subset(
    human_projects,
    project == "SRP009615" & project_type == "data_sources"
)

## Create a RangedSummarizedExperiment (RSE) object at the gene level
rse_gene_SRP009615 <- create_rse(proj_info)

## Explore that RSE object
rse_gene_SRP009615


## Load snapcount R package
library("snapcount")

## snapcount can be used with either a procedural interface
query_jx(compilation = "gtex", regions = "CD99")
query_jx(compilation = "gtex", regions = "CD99", range_filters = samples_count == 10)

## or using the query-builder class
sb <- SnaptronQueryBuilder$new()
sb$compilation("gtex")$regions("CD99")$query_jx()



library("megadepth")

## Next, we locate the example BigWig and annotation files
example_bw <- system.file("tests", "test.bam.all.bw",
    package = "megadepth", mustWork = TRUE
)
annotation_file <- system.file("tests", "testbw2.bed",
    package = "megadepth", mustWork = TRUE
)

## Where are they?
example_bw

## We can then use megadepth to compute the coverage
bw_cov <- get_coverage(
    example_bw,
    op = "mean",
    annotation = annotation_file
)
bw_cov

tissues = c(
"ADIPOSE_TISSUE",
"ADRENAL_GLAND",
"BLADDER",
"BLOOD",
"BLOOD_VESSEL",
"BONE_MARROW",
"BRAIN",
"BREAST",
"CERVIX_UTERI",
"COLON",
"ESOPHAGUS",
"FALLOPIAN_TUBE",
"HEART",
"KIDNEY",
"LIVER",
"LUNG",
"MUSCLE",
"NERVE",
"OVARY",
"PANCREAS",
"PITUITARY",
"PROSTATE",
"SALIVARY_GLAND",
"SKIN",
"SMALL_INTESTINE",
"SPLEEN",
"STOMACH",
"STUDY_NA",
"TESTIS",
"THYROID",
"UTERUS",
"VAGINA")

my_region <- GRanges("chr1", IRanges(100000, 100500))

annotation_file <- "/home/khamanaka/gm_rnaseq/orf/recount/inputs/genes.bed"

for(tissue in tissues){
        gtex_info <- subset(human_projects, file_source == "gtex" & project == tissue & project_type == "data_sources")
        rse_gtex <- create_rse(gtex_info)
        metadata <- colData(rse_gtex)
        bw_urls <- head(rse_gtex$BigWigURL, n=10)
        print(bw_urls)
        res <- get_coverage(bw_urls, op = "mean", annotation = annotation_file)
        print(res)
        # 1. GRanges オブジェクトをデータフレームに変換
        res_df <- as.data.frame(res)
        write.table(res_df, paste0("results/", tissue, ".tsv"), sep = "\t", row.names = F, quote = F)
}
