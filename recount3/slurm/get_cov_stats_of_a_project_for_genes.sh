#!/bin/bash
set -euo pipefail
eval "$(conda shell.bash hook)"
conda activate recount3


PROJECT=$1
GTF=$2
GENE_IDS=$3


Rscript scripts/get_cov_stats_of_a_project_for_genes.R \
        --gtf $GTF \
        --gene_ids $GENE_IDS \
        --project $PROJECT \
        --n_sample 10 \
        --threshold 20 \
        --output results/$PROJECT.txt
