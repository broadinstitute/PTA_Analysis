version 1.0

import "../../tasks/gridss/gridss.wdl" as GridssTools
import "../../tasks/utils/files.wdl" as Utils

workflow get_gridss_vcfs {
    input {
        Array[Pair[String, File]] normal_bams    # [tuple of (sample_id, bam)]
        Array[Pair[String, File]] tumor_bams     # [tuple of (sample_id, bam)]
        String? gridss_driver_vcfs_dir          # Optional directory containing pre-computed GRIDSS VCFs
        String out_dir
    }

    # Combine normal and tumor bams by donor_id
    scatter (normal_bam in normal_bams) {
        String normal_id = normal_bam.left
        scatter (tumor_bam in tumor_bams) {
            if (normal_id == tumor_bam.left) {
                Pair[Pair[String, File], Pair[String, File]] paired_bams = (normal_bam, tumor_bam)
            }
        }
    }
    Array[Pair[Pair[String, File], Pair[String, File]]] input_gridss = select_all(flatten(paired_bams))

    if (defined(gridss_driver_vcfs_dir)) {
        call Utils.extract_gridss_driver_vcf {
            input:
                vcfs_dir = gridss_driver_vcfs_dir
        }

        # Process the extracted VCFs
        scatter (vcf_info in extract_gridss_driver_vcf.vcfs) {
            Array[String] sample_ids = split(vcf_info.sample_id, ";")
            Pair[String, Pair[Array[String], File]] processed_vcf = (
                vcf_info.donor_id,
                (sample_ids, vcf_info.vcf)
            )
        }
    }

    if (!defined(gridss_driver_vcfs_dir)) {
        scatter (bam_pair in input_gridss) {
            call GridssTools.gridss {
                input:
                    normal_bam = bam_pair.left.right,
                    normal_sample_id = bam_pair.left.left,
                    tumor_bam = bam_pair.right.right,
                    tumor_sample_id = bam_pair.right.left,
                    out_dir = out_dir + "/intermediate/svs/gridss/" + bam_pair.left.left
            }
        }
    }

    # Use either pre-computed or newly generated GRIDSS VCFs
    Array[Pair[String, File]] gridss_vcfs = select_first([
        processed_vcf,
        gridss.output_vcf
    ])

    scatter (vcf in gridss_vcfs) {
        call GridssTools.annotate_inserted_sequence {
            input:
                input_vcf = vcf.right,
                donor_id = vcf.left,
                out_dir = out_dir + "/intermediate/svs/gridss/" + vcf.left
        }
    }

    output {
        Array[File] gridss_unfiltered_vcfs = annotate_inserted_sequence.annotated_vcf
    }
}
