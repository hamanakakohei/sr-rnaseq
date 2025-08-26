task rsem_reference {

    File reference_fasta
    File annotation_gtf
    String prefix

    Int memory
    #Int disk_space
    Int num_threads
    #Int num_preempt
    String sge_queue

    command {
        mkdir ${prefix} && cd ${prefix}
        rsem-prepare-reference ${reference_fasta} rsem_reference --gtf ${annotation_gtf} --num-threads ${num_threads}
        cd .. && tar -cvzf ${prefix}.tar.gz ${prefix}
    }

    output {
        File rsem_reference = "${prefix}.tar.gz"
    }

    runtime {
        #docker: "gcr.io/broad-cga-francois-gtex/gtex_rnaseq:V10"
        docker: "gcr.io/broad-cga-francois-gtex/gtex_rnaseq@sha256:80c2db3cec3c08630237665e2d2f044f065022e0bbf7a62d0765f51f811818e2"
        #memory: "${memory}GB"
        memory: "${memory}KB"
        cpu: "${num_threads}"
        #disks: "local-disk ${disk_space} HDD"
        #preemptible: "${num_preempt}"
        sge_queue: "${sge_queue}"
    }

    meta {
        author: "Francois Aguet"
    }
}


workflow rsem_reference_workflow {
    call rsem_reference
}
