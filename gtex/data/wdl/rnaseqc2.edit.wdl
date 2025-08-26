task rnaseqc2 {

    File bam_file
    File genes_gtf
    String sample_id
    String? strandedness
    File? intervals_bed
    File? reference_fasta
    File? reference_fasta_index
    String? flags

    Int memory
    #Int disk_space
    Int num_threads
    #Int num_preempt
    String sge_queue

    command {
        set -euo pipefail
        echo $(date +"[%b %d %H:%M:%S] Running RNA-SeQC 2")
        touch ${sample_id}.fragmentSizes.txt
        touch ${sample_id}.gc_content.tsv
        rnaseqc ${genes_gtf} ${bam_file} . -s ${sample_id} ${"--bed " + intervals_bed} ${"--stranded " + strandedness} ${"--fasta " + reference_fasta} -vv ${flags}
        echo "  * compressing outputs"
        gzip *.gct
        echo $(date +"[%b %d %H:%M:%S] done")
    }

    output {
        File gene_tpm = "${sample_id}.gene_tpm.gct.gz"
        File gene_counts = "${sample_id}.gene_reads.gct.gz"
        File exon_counts = "${sample_id}.exon_reads.gct.gz"
        File metrics = "${sample_id}.metrics.tsv"
        File gc_content = "${sample_id}.gc_content.tsv"
        File insertsize_distr = "${sample_id}.fragmentSizes.txt"
    }

    runtime {
        #docker: "gcr.io/broad-cga-francois-gtex/gtex_rnaseq:V10"
        docker: "gcr.io/broad-cga-francois-gtex/gtex_rnaseq@sha256:80c2db3cec3c08630237665e2d2f044f065022e0bbf7a62d0765f51f811818e2"
        #memory: "${memory}GB"
        memory: "${memory}KB"
        #disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        #preemptible: "${num_preempt}"
        sge_queue: "${sge_queue}"
    }

    #meta {
    #    author: "Francois Aguet"
    #}
}


workflow rnaseqc2_workflow {
    call rnaseqc2
}
