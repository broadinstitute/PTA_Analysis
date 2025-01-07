version 1.0

import "../NextflowModules/Sambamba/0.8.2/Index.wdl" as Sambamba
import "../NextflowModules/GATK/4.2.6.1/GetSampleName.wdl" as GATK

workflow get_indexed_bams {
    input {
        Array[Tuple[String, String, File, File?]] input_bams
    }

    # Filter for unindexed BAMs (those without .bai files)
    Array[Tuple[String, String, File]] unindexed_bams = select_all(
        [for bam in input_bams: 
            if length(bam) == 3 then (bam.left, bam.middle, bam.right) else null]
    )

    # Filter for already indexed BAMs
    Array[Tuple[String, String, File, File]] indexed_bams = select_all(
        [for bam in input_bams: 
            if length(bam) == 4 && basename(bam.right, ".bam") && sub(bam.fourth, ".bai$", "") then bam else null]
    )

    # Index the unindexed BAMs
    scatter (unindexed_bam in unindexed_bams) {
        call Sambamba.Index {
            input:
                bam = unindexed_bam
        }
    }

    # Combine pre-indexed and newly indexed BAMs
    Array[Tuple[String, String, File, File]] all_indexed_bams = flatten([indexed_bams, Index.output_bams])

    # Get sample names for all BAMs
    scatter (indexed_bam in all_indexed_bams) {
        call GATK.GetSampleName {
            input:
                bam = indexed_bam
        }
    }

    output {
        Array[Tuple[String, String, File, File]] output_indexed_bams = GetSampleName.output_bams
    }
}
