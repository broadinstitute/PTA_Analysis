version 1.0

import "../../tasks/bedtools_tasks.wdl" as bedtools
import "../../tasks/htslib_tasks.wdl" as htslib

workflow get_excludeindellist_filtered_vcfs {
    input {
        Array[Pair[String, Pair[String, File]]] somatic_vcfs  # [donor_id, [sample_id, vcf]]
        File? excludeindellist
        String? excludeindellist_filtered_vcfs_dir
        String out_dir
    }

    if (!defined(excludeindellist_filtered_vcfs_dir)) {
        scatter (vcf_pair in somatic_vcfs) {
            String donor_id = vcf_pair.left
            String sample_id = vcf_pair.right.left
            File vcf = vcf_pair.right.right

            call bedtools.intersect as intersectExcludeIndelList {
                input:
                    input_vcf = vcf,
                    exclude_list = select_first([excludeindellist])
            }

            call htslib.bgzip {
                input:
                    input_file = intersectExcludeIndelList.filtered_vcf
            }

            call htslib.tabix {
                input:
                    input_file = bgzip.output_gz
            }
        }
    }

    output {
        Array[Pair[String, Pair[String, Pair[File, File]]]]? processed_vcfs = if !defined(excludeindellist_filtered_vcfs_dir) then 
            zip(
                somatic_vcfs.left,
                zip(
                    somatic_vcfs.right.left,
                    zip(bgzip.output_gz, tabix.output_tbi)
                )
            )
        else none
    }
}