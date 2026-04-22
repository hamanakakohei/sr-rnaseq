#!/bin/bash
#SBATCH --job-name=run_recount3
#SBATCH --partition=epyc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --output=logs/%a.out
#SBATCH --error=logs/%a.err

set -euo pipefail
eval "$(conda shell.bash hook)"
conda activate recount3


PROJECT_LIST=$1 # inputs/projects.list
BED=$2          # inputs/genes.bed
read PROJECT < <(sed -n "${SLURM_ARRAY_TASK_ID}p" "$PROJECT_LIST")


Rscript scripts/get_cov_stats_of_a_project_for_regions.R \
        --annotation $BED \
        --project $PROJECT \
        --n_sample 10 \
        --threshold 20 \
        --output results/$PROJECT.txt \
        #> logs/${SLURM_ARRAY_TASK_ID}.log
