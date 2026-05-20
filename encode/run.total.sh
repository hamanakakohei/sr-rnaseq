# 以下のようにbwをダウンロードして発現量を調べたが、bwが怪しげでQCを色々するはめになったので、BAMの方が楽だったかも
#
# 1. ENCODE portalからマニュアルでデータファイルやメタ情報を取ってくる
# 2. それらのencsrのassemblyやlibrary type、encffのstrandなどを確認する
# 3. それらの情報を1つの表にまとめる
# 4. 表からお気に入りのbwを選んでダウンロード
# 5. bwがstrandが本当に合っているかメジャー遺伝子で確認する
# 6. 各strandのbwのシグナル総数を計算して比が偏り過ぎていないか、合計がどうなっているか（どういう補正をされていたのか）見る
# 7. 興味のあるBED領域に対して補正済みカウント数を計算する
# 8. bwを興味のある領域に絞る
# 9. 全サンプルのbwを足して1つにする

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

awk '$3=="total_RNA-seq"{print $2}' inputs/ENCFF_ENCSR_expType__bwAri.txt \
        | sort -u \
        > inputs/ENCSR__bwAri_total.txt

awk '$3=="total_RNA-seq"{print $1}' inputs/ENCFF_ENCSR_expType__bwAri.txt \
        | sort -u \
        > inputs/ENCFF__bwAri_total.txt


# 2
cat inputs/ENCSR__bwAri_total.txt | while read EXP; do
        scripts/02_check_ENCSR.py --input $EXP
done \
        | sort -u \
        > results/02/encsr_meta__total.txt

cat inputs/ENCFF__bwAri_total.txt | while read EXP; do
        scripts/02_check_ENCFF.py --input $EXP
done \
        | sort -u \
        > results/02/encff_meta__total.txt

awk 'BEGIN{print "encsr expStrand"    }{print $2,     $6}' results/02/encsr_meta__total.txt \
        | sed 's/,//g' \
        > results/02/encsr_expStrand.txt

awk 'BEGIN{print "encff grch bwStrand"}{print $2, $4, $6}' results/02/encff_meta__total.txt \
        | sed 's/,//g' \
        > results/02/encff_bwStrand.txt


# 3
# 注1：以下の2実験はstrand-speだがbwがそうなっていないので除く
# 今回は無し
#
# 注2：以下のencsrはなぜかbwが2種類以上あり、K562とGMというつまらない細胞だったのでもう除いた
# ENCSR000AEE
# ENCSR109IQO
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
#       .query('seqType=="total_RNA-seq" & expStrand!="unstranded" & bwStrand!="signal" & grch=="GRCh38" & encsr!="ENCSR000AEE" & encsr!="ENCSR109IQO"')\
#       .pivot(index="encsr", columns="bwStrand", values="encff")\
#       .reset_index()\
#       .to_csv('results/03/encsr_mins_plus_bw.txt',index=False, sep='\t')


# 4
awk -F"\t" 'NR>1 && $4!="unstranded"{print $7}' results/03/summary.txt \
        | xargs -I{} wget {}
        #| grep -v -e ENCFF743KBR -e ENCFF085IJI \


# 5
# 何サンプルかはigvで見るとchr22がない笑
# 注1：全染色体のデータが揃っているか本当は見た方がよいが、まあ発現高いサンプルだけ見るならいっか、、
# - 今回は無し
#
# 注2：以下はなぜかbwのストランドが逆な気がする --> repaired.txtとして保存した
# - ENCSR571IUZ
# - ENCSR227MWL
# - ENCSR652RSO
# - ENCSR332SRS
# - ENCSR670MXT
# - ENCSR892ZUK
# - ENCSR760INS
# - ENCSR146PLL
# - ENCSR519NFO
# - ENCSR707JVU
# - ENCSR536GUD
# - ENCSR475TVN
# - ENCSR476TKU
# - ENCSR022BYE
# - ENCSR916YUR
# - ENCSR870DRE
# - ENCSR178SNP
# - ENCSR975YGW
# - ENCSR084JIA
# - ENCSR682CFV
# - ENCSR424FAZ
# - ENCSR537BCG
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
                $(awk 'NR>1 && $4!="nan"{sum+=$4}END{print sum}' results/05/${ENCSR}.mins.counts.tab) \
                $(awk 'NR>1 && $5!="nan"{sum+=$5}END{print sum}' results/05/${ENCSR}.mins.counts.tab) \
                $(awk 'NR>1 && $4!="nan"{sum+=$4}END{print sum}' results/05/${ENCSR}.plus.counts.tab) \
                $(awk 'NR>1 && $5!="nan"{sum+=$5}END{print sum}' results/05/${ENCSR}.plus.counts.tab)" \
                >> results/05/stats.txt
done < <(tail -n+2 results/03/encsr_mins_plus_bw.txt)


# 6
# F/R比が一部のサンプルで変
# ENCSR332SRS 0.54462
# ENCSR892ZUK 0.59551
# ENCSR707JVU 0.663615
# ENCSR870IUI 0.670864
# ENCSR000AHH 0.74667
# ENCSR146PLL 0.749698
# ENCSR536GUD 0.749722
# etc.
# ENCSR954PZB 1.59994
# ENCSR434TEU 1.62466
# ENCSR094VRQ 1.64899
# ENCSR158KFO 1.68192
# ENCSR827IXS 1.70094
# ENCSR729CAZ 1.7068
# ENCSR450ENK 1.7228
# ENCSR023ZXN 1.7397
# ENCSR687HJY 1.74945
# ENCSR853WOM 1.7596
# ENCSR146LBD 1.76439
# ENCSR858QEL 1.77432
# ENCSR801MKV 1.80001
# ENCSR504NIU 1.83963
# ENCSR471RUK 1.92807
# ENCSR645TCG 1.9803
# ENCSR802HPM 2.09358
# ENCSR839ZDH 2.14494
# ENCSR257NIR 2.38028
# ENCSR630VJN 2.41196
# ENCSR671WMH 2.42185
# ENCSR436QDU 2.46911
# ENCSR544SAU 2.48646
# ENCSR457ENP 2.54024
# ENCSR150QJY 2.54406
# ENCSR800WIY 2.55853
# ENCSR571RXE 2.68346
# ENCSR080HPT 2.76914
# ENCSR296PMS 2.84378
# ENCSR391VGU 3.8553
# ENCSR403SZN 4.4122
# ENCSR752UNJ 4.43391
# ENCSR450BNZ 6.10878
# ENCSR653ZJF 8.19436
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
done < <(tail -n+430 results/03/encsr_mins_plus_bw.txt)
#done < <(tail -n+2 results/05/encsr_mins_plus_bw.repaired.txt)


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
done < <(tail -n+261 results/03/encsr_mins_plus_bw.txt | awk '{print $2; print $3}')
#done < <(tail -n+2 results/05/encsr_mins_plus_bw.repaired.txt | awk '{print $2; print $3}')


# 8
bwToBg='singularity exec /usr/local/biotools/u/ucsc-bigwigtobedgraph:482--h0b57e2e_0 bigWigToBedGraph'
bgToBw='singularity exec /usr/local/biotools/u/ucsc-bedgraphtobigwig:482--hdc0a859_0 bedGraphToBigWig'
CHROMSIZES=/home/khamanaka/resource/ref/hg38.genome
CHR=chr5
STA=181203000
END=181211000

while read ENCFF; do
        BW=results/04/${ENCFF}.bigWig
        OUT1=results/08/$(basename $BW .bigWig).region.bg
        OUT2=results/08/$(basename $BW .bigWig).region.bw

        $bwToBg \
                -chrom=${CHR} -start=${STA} -end=${END} \
                $BW $OUT1

        $bgToBw \
                $OUT1 $CHROMSIZES $OUT2
done < <(tail -n+2 results/03/encsr_mins_plus_bw.txt | awk '{print $2; print $3}')


# 9
WIGTOOLS='singularity exec /usr/local/biotools/w/wiggletools:1.2.11--h7118728_10 wiggletools'
MERGED_PREFIX=results/09/merged.region

ls -l results/08 | awk '$5==0{print $NF}' | sed 's/.region.bg//' > results/09/no_data_encff.list
awk 'NR>1{print "results/08/"$2".region.bw"}' results/05/encsr_mins_plus_bw.repaired.txt > results/09/bw.f.list
awk 'NR>1{print "results/08/"$3".region.bw"}' results/05/encsr_mins_plus_bw.repaired.txt > results/09/bw.r.list
$WIGTOOLS write_bg ${MERGED_PREFIX}.f.bg sum $(cat results/09/bw.f.list | grep -v -f results/09/no_data_encff.list)
$WIGTOOLS write_bg ${MERGED_PREFIX}.r.bg sum $(cat results/09/bw.r.list | grep -v -f results/09/no_data_encff.list)
$bgToBw ${MERGED_PREFIX}.f.bg $CHROMSIZES ${MERGED_PREFIX}.f.bw
$bgToBw ${MERGED_PREFIX}.r.bg $CHROMSIZES ${MERGED_PREFIX}.r.bw
