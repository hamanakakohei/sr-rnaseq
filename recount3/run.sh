#!/usr/bin/env bash
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
