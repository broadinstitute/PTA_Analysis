version 1.0

import "../../tasks/randomForest.wdl" as RF
import "../../tasks/htslib.wdl" as htslib
import "../../tasks/utils.wdl" as utils

workflow get_ptato_vcfs {
    input {
        Array[Pair[String, Pair[String, File]]] somatic_vcfs  # [donor_id, [sample_id, vcf]]
        Array[Pair[String, Pair[String, File]]] rf_tables     # [donor_id, [sample_id, rf_table]]
        String? ptato_vcfs_dir
        String out_dir
    }

    if (defined(ptato_vcfs_dir)) {
        call utils.extractPtatoVcfFromDir {
            input:
                directory = select_first([ptato_vcfs_dir])
        }
        call htslib.bgzip as gzip_extracted {
            input:
                input_file = extractPtatoVcfFromDir.vcf
        }
    }

    if (!defined(ptato_vcfs_dir)) {
        scatter (vcf_pair in somatic_vcfs) {
            String donor_id = vcf_pair.left
            String sample_id = vcf_pair.right.left
            File vcf = vcf_pair.right.right

            # Find matching RF table
            File rf_table = select_first([rf_tables[donor_id][sample_id]])

            call RF.test_snv_rf {
                input:
                    vcf = vcf,
                    rf_table = rf_table
            }

            call htslib.bgzip {
                input:
                    input_file = test_snv_rf.filtered_vcf
            }

            call htslib.tabix {
                input:
                    vcf_gz = bgzip.output_gz
            }

            # Copy files to output directory
            String vcf_name = basename(bgzip.output_gz)
            String tbi_name = basename(tabix.output_tbi)
            String output_path = out_dir + "/snvs/" + donor_id

            call utils.copy_files {
                input:
                    source_vcf = bgzip.output_gz,
                    source_tbi = tabix.output_tbi,
                    destination_dir = output_path,
                    vcf_filename = vcf_name,
                    tbi_filename = tbi_name
            }
        }
    }

    output {
        Array[Pair[String, Pair[String, Pair[File, File]]]] ptato_vcfs = 
            if defined(ptato_vcfs_dir) then select_all([extractPtatoVcfFromDir.vcf_pairs])
            else select_all(copy_files.output_pairs)
    }
}
