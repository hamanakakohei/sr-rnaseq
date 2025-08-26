task markduplicates {

    File edited_script
    File input_bam
    String prefix
    Int? max_records_in_ram
    Float? sorting_collection_size_ratio

    Float memory
    Int num_threads
    #Int java_memory = floor(memory - 0.5)
    Int java_memory = floor(memory - 0.5) * num_threads
    #Int disk_space
    #Int num_preempt
    Int ParallelGCThreads
    String sge_queue

    String output_bam = sub(basename(input_bam), "\\.bam$", ".md.bam")

    command {
        set -euo pipefail
        #python3 -u /src/run_MarkDuplicates.py ${input_bam} ${prefix} \
        python3 -u ${edited_script} ${input_bam} ${prefix} \
            --ParallelGCThreads ${ParallelGCThreads} \
            --memory ${java_memory} \
            ${"--max_records_in_ram " + max_records_in_ram} \
            ${"--sorting_collection_size_ratio " + sorting_collection_size_ratio}
        samtools index ${output_bam}
    }

    output {
        File bam_file = "${output_bam}"
        File bam_index = "${output_bam}.bai"
        File metrics = "${prefix}.marked_dup_metrics.txt"
    }

    runtime {
        #docker: "gcr.io/broad-cga-francois-gtex/gtex_rnaseq:V10"
        docker: "gcr.io/broad-cga-francois-gtex/gtex_rnaseq@sha256:80c2db3cec3c08630237665e2d2f044f065022e0bbf7a62d0765f51f811818e2"
        #memory: "${memory}GB"
        memory: "${memory}KB"
        #disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        #preemptible: "${num_preempt}"
        ParallelGCThreads: "${ParallelGCThreads}"
        sge_queue: "${sge_queue}"
    }

    meta {
        author: "Francois Aguet"
    }
}


workflow markduplicates_workflow {
    call markduplicates
}
