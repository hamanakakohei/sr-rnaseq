#import "https://api.firecloud.org/ga4gh/v1/tools/broadinstitute_gtex:fastqc_v1-0_BETA/versions/2/plain-WDL/descriptor" as fastqc_wdl
#import "https://api.firecloud.org/ga4gh/v1/tools/broadinstitute_gtex:star_v1-0_BETA/versions/8/plain-WDL/descriptor" as star_wdl
#import "https://api.firecloud.org/ga4gh/v1/tools/broadinstitute_gtex:markduplicates_v1-0_BETA/versions/6/plain-WDL/descriptor" as markduplicates_wdl
#import "https://api.firecloud.org/ga4gh/v1/tools/broadinstitute_gtex:rsem_v1-0_BETA/versions/6/plain-WDL/descriptor" as rsem_wdl
#import "https://api.firecloud.org/ga4gh/v1/tools/broadinstitute_gtex:rnaseqc2_v1-0_BETA/versions/4/plain-WDL/descriptor" as rnaseqc_wdl
import "/path/to/star.edit.wdl" as star_wdl
import "/path/to/markduplicates.edit.wdl" as markduplicates_wdl
import "/path/to/rsem.edit.wdl" as rsem_wdl
import "/path/to/rnaseqc2.edit.wdl" as rnaseqc_wdl

workflow rnaseq_pipeline_fastq_workflow {

    String prefix

    #call fastqc_wdl.fastqc {}

    call star_wdl.star {
        input: prefix=prefix
    }

    call markduplicates_wdl.markduplicates {
        input: input_bam=star.bam_file, prefix=prefix
    }

    call rsem_wdl.rsem {
        input: transcriptome_bam=star.transcriptome_bam, prefix=prefix
    }

    call rnaseqc_wdl.rnaseqc2 {
        input: bam_file=markduplicates.bam_file, sample_id=prefix
    }
}
