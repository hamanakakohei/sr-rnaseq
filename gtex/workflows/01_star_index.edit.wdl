task star_index {

    File reference_fasta
    File annotation_gtf
    String prefix
    Int overhang
    Int? suffix_length_max
    String? transform_type
    File? transform_vcf

    Int memory
    #Int disk_space
    Int num_threads
    #Int num_preempt
    String sge_queue

    command {
        set -euo pipefail
        mkdir ${prefix}
        STAR \
            --runMode genomeGenerate \
            --genomeDir ${prefix} \
            --genomeFastaFiles ${reference_fasta} \
            --sjdbGTFfile ${annotation_gtf} \
            --sjdbOverhang ${overhang} \
            ${"--genomeSuffixLengthMax " + suffix_length_max} \
            ${"--genomeTransformType " + transform_type} \
            ${"--genomeTransformVCF " + transform_vcf} \
            --runThreadN ${num_threads}
        tar -cvzf ${prefix}.tar.gz ${prefix}
    }

    output {
        File star_index = "${prefix}.tar.gz"
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

    #meta {
    #    author: "Francois Aguet"
    #}
}


workflow star_index_workflow {
    call star_index
}
