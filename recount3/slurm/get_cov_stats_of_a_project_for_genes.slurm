#!/bin/bash
#SBATCH --job-name=run_recount3
#SBATCH --partition=medium
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --output=logs/%a.out
#SBATCH --error=logs/%a.err

set -euo pipefail
eval "$(conda shell.bash hook)"
conda activate recount3


PROJECT_LIST=$1  # inputs/projects.list
GTF=$2
GENE_IDS=$3
read PROJECT < <(sed -n "${SLURM_ARRAY_TASK_ID}p" "$PROJECT_LIST")


Rscript scripts/get_cov_stats_of_a_project_for_genes.R \
        --gtf $GTF \
        --gene_ids $GENE_IDS \
        --project $PROJECT \
        --n_sample 10 \
        --threshold 20 \
        --output results/$PROJECT.txt
