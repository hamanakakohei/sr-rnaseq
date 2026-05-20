#!/usr/bin/env bash
# 
# 1：各研究の各実験（数が多すぎると適当な数に絞る）について、指定した遺伝子領域のシグナルを補正しつつ合計する
# 2：ホームページからダウンロードしたメタ情報（研究レベル）を1の結果と結合する
# 3：gtexとtcgaについては組織毎にシグナルを平均する
set -euo pipefail


# 1
PROJECT_LIST=inputs/projects.list
#BED=inputs/genes.bed
GTF=inputs/gencodev47.gtf
GENE_IDS=inputs/gene_ids.list


# slurm ver.
STA=1
END=$(wc -l < $PROJECT_LIST)
#END=1000

sbatch \
         --array=${STA}-${END}%50 \
         slurm/get_cov_stats_of_a_project_for_genes.slurm \
         $PROJECT_LIST \
         $GTF \
         $GENE_IDS
         #--array=1,3,7,8,9%10 \
         #--array=50%10 \


## local ver.
#while read PROJECT; do
#       echo $PROJECT
#       ./slurm/get_cov_stats_of_a_project_for_genes.sh \
#               $PROJECT \
#               $GTF \
#               $GENE_IDS \
#               > logs/local/$PROJECT.log 2>&1
#done < <(tail -n+1 $PROJECT_LIST)

# local ver. 2
tail -n+1    $PROJECT_LIST | head -n 3000 | tail -n+550 | xargs -P 20 -I {} ./slurm/get_cov_stats_of_a_project_for_genes.sh {} $GTF $GENE_IDS
tail -n+3001 $PROJECT_LIST | head -n 3000 |               xargs -P 20 -I {} ./slurm/get_cov_stats_of_a_project_for_genes.sh {} $GTF $GENE_IDS
tail -n+6001 $PROJECT_LIST                |               xargs -P 20 -I {} ./slurm/get_cov_stats_of_a_project_for_genes.sh {} $GTF $GENE_IDS


# 2
META=inputs/recount3_selection_2026-05-11_13_55_57.215489.csv

./scripts/merge_results.py \
    --results_dir results \
    --csv $META \
    --output merged_results.tsv



