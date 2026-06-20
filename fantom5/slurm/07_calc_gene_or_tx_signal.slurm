#!/bin/bash
#SBATCH --job-name=07
#SBATCH --partition=rome
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --output=07_10_%a.out
#SBATCH --error=07_10_%a.err
set -euo pipefail
eval "$(conda shell.bash hook)"
conda activate misc_20250301


GTF=inputs/gtf
IDS_OI=inputs/gene_ids_oi.list
TSS_MARGIN=10

mkdir -p results_${TSS_MARGIN}/07/

while read ENCFF; do
        echo "Checking ${ENCFF}..."
        BW=inputs/07/${ENCFF}.bw

        python3 scripts/07_calc_gene_or_tx_signal.py \
                --bw $BW \
                --gtf $GTF \
                --mode gene \
                --list $IDS_OI \
                --tss \
                --tss_interval ${TSS_MARGIN},${TSS_MARGIN} \
                --out results_${TSS_MARGIN}/07/${ENCFF}.TotalSignal.txt
done < <(cat inputs/07/sample_fwd_rev_bw.list | awk '{print $2; print $3}')
