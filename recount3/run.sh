#!/usr/bin/env bash
set -euo pipefail


PROJECT_LIST=inputs/projects.list
BED=inputs/genes.bed


STA=1
END=$(wc -l < $PROJECT_LIST)
#END=1000


sbatch \
         --array=${STA}-${END}%1 \
         slurm/get_cov_stats_of_a_project_for_regions.slurm \
         $PROJECT_LIST \
         $BED
         #--array=1,3,7,8,9%10 \
         #--array=50%10 \
