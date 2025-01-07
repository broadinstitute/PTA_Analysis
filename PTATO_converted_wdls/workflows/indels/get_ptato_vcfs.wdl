version 1.0

# Import tasks from other WDL files
import "../../tasks/Utils/excludeIndels.wdl" as ExcludeIndels
import "../../tasks/htslib/bgzip.wdl" as Bgzip
import "../../tasks/htslib/tabix.wdl" as Tabix
import "../../tasks/Utils/extractPtatoVcf.wdl" as ExtractPtato

workflow get_ptato_vcfs {
    input {
        Array[Pair[String, Pair[String, File]]] somatic_vcfs  # [donor_id, [sample_id, vcf]]
        Array[Pair[String, Pair[String, File]]] context_beds  # [donor_id, [sample_id, bed]]
        String? ptato_vcfs_dir
        File exclude_indel_list
        String out_dir
    }

    if (defined(ptato_vcfs_dir)) {
        call ExtractPtato.extractPtatoVcfFromDir {
            input:
                directory = select_first([ptato_vcfs_dir])
        }
    }

    if (!defined(ptato_vcfs_dir)) {
        scatter (vcf_pair in somatic_vcfs) {
            String donor_id = vcf_pair.left
            Pair[String, File] sample_vcf = vcf_pair.right
            
            # Find matching context bed
            Pair[String, File] context_bed = select_first([
                select_all(context_beds.right[context_beds.left == donor_id])
            ])

            call ExcludeIndels.excludeIndels {
                input:
                    vcf = sample_vcf.right,
                    context_bed = context_bed.right,
                    exclude_list = exclude_indel_list
            }

            call Bgzip.bgzip {
                input:
                    input_file = excludeIndels.filtered_vcf
            }

            call Tabix.tabix {
                input:
                    vcf_gz = bgzip.compressed_file
            }

            # Copy files to output directory
            String vcf_name = basename(bgzip.compressed_file)
            String tbi_name = basename(tabix.index_file)
            String donor_out_dir = out_dir + "/indels/" + donor_id

            call copy_files {
                input:
                    vcf_gz = bgzip.compressed_file,
                    vcf_tbi = tabix.index_file,
                    destination_dir = donor_out_dir,
                    vcf_name = vcf_name,
                    tbi_name = tbi_name
            }
        }
    }

    output {
        Array[Pair[File, File]] ptato_vcfs = select_first([
            extractPtatoVcfFromDir.ptato_vcfs,
            copy_files.copied_files
        ])
    }
}

task copy_files {
    input {
        File vcf_gz
        File vcf_tbi
        String destination_dir
        String vcf_name
        String tbi_name
    }

    command <<<
        mkdir -p ~{destination_dir}
        cp ~{vcf_gz} ~{destination_dir}/~{vcf_name}
        cp ~{vcf_tbi} ~{destination_dir}/~{tbi_name}
    >>>

    output {
        Pair[File, File] copied_files = (
            destination_dir + "/" + vcf_name,
            destination_dir + "/" + tbi_name
        )
    }

    runtime {
        docker: "ubuntu:latest"
    }
}