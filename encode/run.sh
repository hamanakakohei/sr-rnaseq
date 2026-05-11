# 以下のようにbwをダウンロードして発現量を調べたが、bwが怪しげでQCを色々するはめになったので、BAMの方が楽だったかも
#
# 1. ENCODE portalからマニュアルでデータファイルやメタ情報を取ってくる
# 2. それらのencsrのassemblyやlibrary type、encffのstrandなどを確認する
# 3. それらの情報を1つの表にまとめる
# 4. 表からお気に入りのbwを選んでダウンロード
# 5. bwがstrandが本当に合っているかメジャー遺伝子で確認する
# 6. 各strandのbwのシグナル総数を計算して比が偏り過ぎていないか、合計がどうなっているか（どういう補正をされていたのか）見る
# 7. 興味のあるBED領域に対して補正済みカウント数を計算する

# 1
cut -f2,4,16,60,75,76,77 inputs/experiment_report_2026_4_22_5h_14m.tsv \
        | tail -n+3 \
        | tr ' ' '_' \
        | awk -F"[\t,]" '{for(i=3;i<=NF;i++) print $i"\t"$1"\t"$2}' \
        | awk -F"\t" '$1!=""' \
        | sed 's/files//g; s/\///g' \
        | sort -u \
        > inputs/ENCFF_ENCSR_expType.txt

awk -F"/" 'NR>1{print $NF}' inputs/files_default_files.txt \
        | sed 's/.bigWig//' \
        | sort -u \
        > inputs/ENCFF.txt

join inputs/ENCFF.txt inputs/ENCFF_ENCSR_expType.txt \
        > inputs/ENCFF_ENCSR_expType__bwAri.txt

awk '$3=="polyA_plus_RNA-seq"{print $2}' inputs/ENCFF_ENCSR_expType__bwAri.txt \
        > inputs/ENCSR__bwAri_polyA.txt

awk '$3=="polyA_plus_RNA-seq"{print $1}' inputs/ENCFF_ENCSR_expType__bwAri.txt \
        > inputs/ENCFF__bwAri_polyA.txt

awk '$3=="total_RNA-seq"{print $2}' inputs/ENCFF_ENCSR_expType__bwAri.txt \
        > inputs/ENCSR__bwAri_total.txt

awk '$3=="total_RNA-seq"{print $1}' inputs/ENCFF_ENCSR_expType__bwAri.txt \
        > inputs/ENCFF__bwAri_total.txt


# 2
cat inputs/ENCSR__bwAri_polyA.txt | while read EXP; do
        scripts/02_check_ENCSR.py --input $EXP
done > results/02/encsr_meta__polyA.txt

cat inputs/ENCFF__bwAri_polyA.txt | while read EXP; do
        scripts/02_check_ENCFF.py --input $EXP
done > results/02/encff_meta__polyA.txt

cat inputs/ENCSR__bwAri_total.txt | while read EXP; do
        scripts/02_check_ENCSR.py --input $EXP
done > results/02/encsr_meta__total.txt

cat inputs/ENCFF__bwAri_total.txt | while read EXP; do
        scripts/02_check_ENCFF.py --input $EXP
done > results/02/encff_meta__total.txt

awk 'BEGIN{print "encsr expStrand"    }{print $2,     $6}' results/02/encsr_meta__polyA.txt \
        | sed 's/,//g' \
        > results/02/encsr_expStrand.txt

awk 'BEGIN{print "encff grch bwStrand"}{print $2, $4, $6}' results/02/encff_meta__polyA.txt \
        | sed 's/,//g' \
        > results/02/encff_bwStrand.txt


# 3
# 注1：以下の2実験はstrand-speだがbwがそうなっていないので除く
# bwStrand  strand       size
# minus     forward        52
#           reverse       135
# plus      forward        52
#           reverse       135
# signal    forward         1
#           reverse         1
#           unstranded    208
#
#       encsr        encff bwStrand   strand
# ENCSR000BZV  ENCFF743KBR   signal  forward
# ENCSR783BUO  ENCFF085IJI   signal  reverse
#
# 注2：以下のencsrはなぜかbwが2種類あり片方を除いた
# ENCSR282KJZ ENCFF145MJJ plus: 431Mb
# ENCSR282KJZ ENCFF369YHB plus: 405Mb
# ENCSR282KJZ ENCFF981SRF minus: 405Mb
# ENCSR282KJZ ENCFF962XNB minus: 380Mb
#
#
# import pandas as pd
#
# encsr_expStrand = pd.read_table('results/02/encsr_expStrand.txt', sep=' ')
# encff_bwStrand  = pd.read_table('results/02/encff_bwStrand.txt',  sep=' ')
# encff_encsr_seqType = pd.read_table("inputs/ENCFF_ENCSR_expType__bwAri.txt", sep=" ", names=['encff','encsr','seqType'])
# url_encff = pd.read_table('inputs/files_default_files.txt', skiprows=1, names=['url'])
# url_encff['encff'] = url_encff['url'].map(lambda x: x.split('/')[4])
# encsr_meta = pd.read_table("inputs/experiment_report_2026_4_22_5h_14m.tsv", skiprows=1).\
#       rename({'Accession': 'encsr'}, axis=1)
#
# pd.merge(
#       pd.merge(
#               pd.merge(
#                       pd.merge(encff_encsr_seqType, encsr_expStrand),
#                       encff_bwStrand),
#               url_encff),
#       encsr_meta
#       )\
#       .drop_duplicates()\
#       .to_csv('results/03/summary.txt', sep='\t', index=False)
#
# pd.read_table('results/03/summary.txt')\
#       .query('seqType=="polyA_plus_RNA-seq" & expStrand!="unstranded" & bwStrand!="signal" & encff!="ENCFF369YHB" & encff!="ENCFF962XNB"')\
#       .pivot(index="encsr", columns="bwStrand", values="encff")\
#       .reset_index()\
#       .to_csv('results/03/encsr_mins_plus_bw.txt',index=False, sep='\t')


# 4
awk -F"\t" 'NR>1 && $4!="unstranded"{print $7}' results/03/summary.txt \
        | grep -v -e ENCFF743KBR -e ENCFF085IJI \
        | xargs -I{} wget {}


# 5
# 何サンプルかはigvで見るとchr22がない笑
# 注1：全染色体のデータが揃っているか本当は見た方がよいが、まあ発現高いサンプルだけ見るならいっか、、
# - ENCSR196ARY
# - ENCSR329ZRF
# - ENCSR445DAC
# - ENCSR920UAO
#
# 注2：以下はなぜかbwのストランドが逆な気がする --> repaired.txtとして保存した
# - ENCSR000BZU
# - ENCSR270XRV
# - ENCSR552EGO
# - ENCSR329MHM
# - ENCSR643QIZ
# - ENCSR297UBP
# - ENCSR654UPQ
# - ENCSR637VLS
GENCODE_PLUS_BED=/home/khamanaka/resource/gencode/gencode.v47.annotation.plus.chr22.bed
GENCODE_MINS_BED=/home/khamanaka/resource/gencode/gencode.v47.annotation.mins.chr22.bed

while read ENCSR MINS PLUS; do
        echo "Checking ${ENCSR}..."
        PLUS_BW=results/04/${PLUS}.bigWig
        MINS_BW=results/04/${MINS}.bigWig

        # プラス鎖遺伝子領域でのプラス/マイナスBigWigの平均値を計算
        multiBigwigSummary BED-file \
                -b ${PLUS_BW} ${MINS_BW} \
                --BED ${GENCODE_PLUS_BED} \
                -out results/05/${ENCSR}.plus.npz \
                --outRawCounts results/05/${ENCSR}.plus.counts.tab

        multiBigwigSummary BED-file \
                -b ${PLUS_BW} ${MINS_BW} \
                --BED ${GENCODE_MINS_BED} \
                -out results/05/${ENCSR}.mins.npz \
                --outRawCounts results/05/${ENCSR}.mins.counts.tab

        echo -e \
                "${ENCSR} \
                $(awk 'NR>1 && $4!="nan"{sum+=$4}END{print sum}' results/04/${ENCSR}.mins.counts.tab) \
                $(awk 'NR>1 && $5!="nan"{sum+=$5}END{print sum}' results/04/${ENCSR}.mins.counts.tab) \
                $(awk 'NR>1 && $4!="nan"{sum+=$4}END{print sum}' results/04/${ENCSR}.plus.counts.tab) \
                $(awk 'NR>1 && $5!="nan"{sum+=$5}END{print sum}' results/04/${ENCSR}.plus.counts.tab)" \
                >> results/05/stats.txt
done < <(tail -n+2 results/03/encsr_mins_plus_bw.txt)


# 6
# F/R比が一部のサンプルで変
# ENCSR701UNO 0.640044
# ENCSR920UAO 0.72274
# ENCSR329ZRF 0.729767
# ENCSR249CKG 0.819849
# ENCSR192NBO 0.847022
# etc.
# ENCSR675YAS 1.63684
# ENCSR518XRJ 1.63975
# ENCSR085HNI 1.66571
# ENCSR670WQY 1.69224
# ENCSR775KCE 1.69456
# ENCSR612HYR 1.69895
# ENCSR910QOX 1.72529
# ENCSR995BHD 1.7741
# ENCSR001UXR 1.84701
# ENCSR510PSL 1.92136
# ENCSR571BML 2.0158
# ENCSR071ZMO 2.15115
# ENCSR629VMZ 2.15807
# ENCSR922VBO 2.29563
# ENCSR769LNJ 2.41952
# ENCSR635GTY 2.57406
# ENCSR146ZKR 2.59216
# ENCSR270OKS 3.00839
# ENCSR433XCV 3.92116
# ENCSR598KJX 4.51559
# ENCSR721HDG 4.61297
# ENCSR741QEH 6.31748
# ENCSR763NOO 7.058
# ENCSR719HRO 8.67672
# ENCSR980UEY 9.16283
# ENCSR039ICU 10.6666
# ENCSR999ZCI 11.1459
# ENCSR825GWD 12.8154
APP=/usr/local/biotools/u/ucsc-bigwigtobedgraph\:469--h664eb37_1
while read ENCSR MINS PLUS; do
        echo "Checking ${ENCSR}..."
        PLUS_BW=results/04/${PLUS}.bigWig
        MINS_BW=results/04/${MINS}.bigWig
        PLUS_BG=results/06/tmp.plus.bg
        MINS_BG=results/06/tmp.mins.bg

        singularity exec $APP \
                bigWigToBedGraph \
                $PLUS_BW \
                $PLUS_BG

        singularity exec $APP \
                bigWigToBedGraph \
                $MINS_BW \
                $MINS_BG

        echo -e \
                "${ENCSR} \
                $(awk '{sum += ($3 - $2) * $4} END {print sum}' $PLUS_BG ) \
                $(awk '{sum += ($3 - $2) * $4} END {print sum}' $MINS_BG )" \
                >> results/06/stats.txt
done < <(tail -n+176 results/05/encsr_mins_plus_bw.repaired.txt)
#done < <(tail -n+176 results/03/encsr_mins_plus_bw.txt)


# 7
# 以下の2サンプルは恐らく特定の染色体領域ノデータがない？でエラーになる
# - ENCSR445DAC
# - ENCSR920UA
GTF=inputs/07/gm24_gencodev47.chr.scaffold.444141isoforms.geneid_corrected.geneid_refined.add_txType_geneType_geneRow.sort.gtf
IDS_OI=inputs/07/gene_ids_oi.list

while read ENCFF; do
        echo "Checking ${ENCFF}..."
        BW=results/04/${ENCFF}.bigWig

        scripts/07_calc_gene_or_tx_signal.py \
                --bw $BW \
                --gtf $GTF \
                --mode gene \
                --list $IDS_OI \
                --out results/07/${ENCFF}.TotalSignal.txt
done < <(tail -n+155 results/05/encsr_mins_plus_bw.repaired.txt | awk '{print $2; print $3}')


# 8
scripts/build_expression_matrix.py \
        --expr-dir results/07/ \
        --mapping results/05/encsr_mins_plus_bw.repaired.txt \
        --summary results/03/summary.txt \
        --output results/08/expression_matrix.tsv
