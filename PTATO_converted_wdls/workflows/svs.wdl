version 1.0

# Import required tasks
import "../tasks/get_gzipped_vcfs.wdl" as gzipped_vcfs
import "../tasks/svs/get_gridss_vcfs.wdl" as gridss
import "../tasks/svs/get_gripss_vcfs.wdl" as gripss
import "../tasks/svs/combine_cnvs_svs.wdl" as combine
import "../tasks/utils/get_files_from_dir.wdl" as utils

workflow svs {
    input {
        Array[Pair[String, Pair[String, File]]] bams          # [donor_id, [sample_id, bam]]
        Array[Pair[String, Pair[String, File]]] germline_vcfs # [donor_id, [sample_id, vcf]]
        Array[Pair[String, Pair[String, File]]] filtered_cnv_files
        Array[String] bulk_names
        String? gripss_somatic_filtered_vcfs_dir
        String? gridss_unfiltered_vcfs_dir
    }

    # Separate normal and tumor BAMs
    scatter (bam in bams) {
        String donor = bam.left
        String sample = bam.right.left
        File bam_file = bam.right.right

        if (sample == "Normal") {
            Pair[String, Pair[String, File]] normal_bam = bam
        }
        if (sample != "Normal") {
            Pair[String, Pair[String, File]] tumor_bam = bam
        }
    }

    Array[Pair[String, Pair[String, File]]] normal_bams = select_all(normal_bam)
    Array[Pair[String, Pair[String, File]]] tumor_bams = select_all(tumor_bam)

    if (defined(gripss_somatic_filtered_vcfs_dir)) {
        call utils.extractGripssSomaticFilteredVcfFromDir {
            input:
                directory = select_first([gripss_somatic_filtered_vcfs_dir])
        }
        
        call gzipped_vcfs.get_gzipped_vcfs as gripss_vcfs {
            input:
                vcfs = extractGripssSomaticFilteredVcfFromDir.vcfs
        }
    }

    if (!defined(gripss_somatic_filtered_vcfs_dir)) {
        if (defined(gridss_unfiltered_vcfs_dir)) {
            call utils.extractGridssUnfilteredVcfFromDir {
                input:
                    directory = select_first([gridss_unfiltered_vcfs_dir])
            }
            
            call gzipped_vcfs.get_gzipped_vcfs as gridss_raw_vcfs {
                input:
                    vcfs = extractGridssUnfilteredVcfFromDir.vcfs
            }
        }

        if (!defined(gridss_unfiltered_vcfs_dir)) {
            call gridss.get_gridss_vcfs {
                input:
                    normal_bams = normal_bams,
                    tumor_bams = tumor_bams
            }
        }

        File gridss_vcfs = select_first([gridss_raw_vcfs.output_vcfs, get_gridss_vcfs.output_vcfs])

        call gripss.get_gripss_vcfs {
            input:
                gridss_vcfs = gridss_vcfs
        }
    }

    File final_gripss_vcfs = select_first([gripss_vcfs.output_vcfs, get_gripss_vcfs.output_vcfs])

    call combine.combine_cnvs_svs {
        input:
            cnv_files = filtered_cnv_files,
            sv_files = final_gripss_vcfs
    }

    output {
        File combined_variants = combine_cnvs_svs.output_file
    }
}
