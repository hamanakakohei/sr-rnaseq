task rsem {

    File transcriptome_bam
    File rsem_reference
    String prefix

    Int memory
    #Int disk_space
    Int num_threads
    #Int num_preempt

    Int? max_frag_len
    String? estimate_rspd
    String? is_stranded
    String? paired_end
    String sge_queue

    command {
        set -euo pipefail
        mkdir rsem_reference
        tar -xvvf ${rsem_reference} -C rsem_reference --strip-components=1

        /src/run_RSEM.py \
            ${"--max_frag_len " + max_frag_len} \
            ${"--estimate_rspd " + estimate_rspd} \
            ${"--is_stranded " + is_stranded} \
            ${"--paired_end " + paired_end} \
            --threads ${num_threads} \
            rsem_reference ${transcriptome_bam} ${prefix}
        gzip *.results
    }

    output {
        File genes="${prefix}.rsem.genes.results.gz"
        File isoforms="${prefix}.rsem.isoforms.results.gz"
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


workflow rsem_workflow {
    call rsem
}
