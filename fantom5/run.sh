#!/usr/bin/en bash
# encode repoのスクリプトを流用する
#
# サンプル数：
#     5 human.cell_line.LQhCAGE
#   266 human.cell_line.hCAGE
#    21 human.fractionation.hCAGE
#    52 human.primary_cell.LQhCAGE
#   512 human.primary_cell.hCAGE
#     2 human.timecourse.LQhCAGE
#   783 human.timecourse.hCAGE
#   188 human.tissue.hCAGE


# 7
GTF=inputs/gtf
IDS_OI=inputs/gene_ids_oi.list
TSS_MARGIN=100

sbatch \
        slurm/07_calc_gene_or_tx_signal.slurm \
        $GTF \
        $IDS_OI \
        $TSS_MARGIN


# 8
ls human.cell_line.hCAGE/results/07      | sed -e 's/.fwd.TotalSignal.txt//; s/.rev.TotalSignal.txt//' | sort -u | awk 'BEGIN{print "encsr"}{print}' > human.cell_line.hCAGE/results/03/summary.txt
ls human.cell_line.LQhCAGE/results/07    | sed -e 's/.fwd.TotalSignal.txt//; s/.rev.TotalSignal.txt//' | sort -u | awk 'BEGIN{print "encsr"}{print}' > human.cell_line.LQhCAGE/results/03/summary.txt
ls human.cell_line.hCAGE/results/07      | sed -e 's/.fwd.TotalSignal.txt//; s/.rev.TotalSignal.txt//' | sort -u | awk 'BEGIN{print "encsr"}{print}' > human.cell_line.hCAGE/results/03/summary.txt
ls human.fractionation.hCAGE/results/07  | sed -e 's/.fwd.TotalSignal.txt//; s/.rev.TotalSignal.txt//' | sort -u | awk 'BEGIN{print "encsr"}{print}' > human.fractionation.hCAGE/results/03/summary.txt
ls human.primary_cell.LQhCAGE/results/07 | sed -e 's/.fwd.TotalSignal.txt//; s/.rev.TotalSignal.txt//' | sort -u | awk 'BEGIN{print "encsr"}{print}' > human.primary_cell.LQhCAGE/results/03/summary.txt
ls human.primary_cell.hCAGE/results/07   | sed -e 's/.fwd.TotalSignal.txt//; s/.rev.TotalSignal.txt//' | sort -u | awk 'BEGIN{print "encsr"}{print}' > human.primary_cell.hCAGE/results/03/summary.txt
ls human.timecourse.hCAGE/results/07     | sed -e 's/.fwd.TotalSignal.txt//; s/.rev.TotalSignal.txt//' | sort -u | awk 'BEGIN{print "encsr"}{print}' > human.timecourse.hCAGE/results/03/summary.txt
ls human.tissue.hCAGE/results/07         | sed -e 's/.fwd.TotalSignal.txt//; s/.rev.TotalSignal.txt//' | sort -u | awk 'BEGIN{print "encsr"}{print}' > human.tissue.hCAGE/results/03/summary.txt

awk 'NR==1{print "encsr\tminus\tplus"}NR>1{print $1"\t"$1".rev\t"$1".fwd"}' human.cell_line.LQhCAGE/results/03/summary.txt    > human.cell_line.LQhCAGE/results/05/encsr_mins_plus_bw.repaired.txt
awk 'NR==1{print "encsr\tminus\tplus"}NR>1{print $1"\t"$1".rev\t"$1".fwd"}' human.cell_line.hCAGE/results/03/summary.txt      > human.cell_line.hCAGE/results/05/encsr_mins_plus_bw.repaired.txt
awk 'NR==1{print "encsr\tminus\tplus"}NR>1{print $1"\t"$1".rev\t"$1".fwd"}' human.fractionation.hCAGE/results/03/summary.txt  > human.fractionation.hCAGE/results/05/encsr_mins_plus_bw.repaired.txt
awk 'NR==1{print "encsr\tminus\tplus"}NR>1{print $1"\t"$1".rev\t"$1".fwd"}' human.primary_cell.LQhCAGE/results/03/summary.txt > human.primary_cell.LQhCAGE/results/05/encsr_mins_plus_bw.repaired.txt
awk 'NR==1{print "encsr\tminus\tplus"}NR>1{print $1"\t"$1".rev\t"$1".fwd"}' human.primary_cell.hCAGE/results/03/summary.txt   > human.primary_cell.hCAGE/results/05/encsr_mins_plus_bw.repaired.txt
awk 'NR==1{print "encsr\tminus\tplus"}NR>1{print $1"\t"$1".rev\t"$1".fwd"}' human.timecourse.hCAGE/results/03/summary.txt     > human.timecourse.hCAGE/results/05/encsr_mins_plus_bw.repaired.txt
awk 'NR==1{print "encsr\tminus\tplus"}NR>1{print $1"\t"$1".rev\t"$1".fwd"}' human.tissue.hCAGE/results/03/summary.txt         > human.tissue.hCAGE/results/05/encsr_mins_plus_bw.repaired.txt

scripts/08_build_expresion_matrix.py \
        --expr-dir results/07/ \
        --mapping results/05/encsr_mins_plus_bw.repaired.txt \
        --summary results/03/summary.txt \
        --output results/08/expression_matrix.tsv
