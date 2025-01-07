version 1.0

# Import tasks from other WDL files
import "../../tasks/ptato_filter.wdl" as PtatoFilter
import "../../tasks/htslib.wdl" as Htslib

workflow filter_ptato_vcfs {
    input {
        Array[Pair[Pair[String, String], File]] ptato_vcfs  # Array of [donor_id, sample_id, vcf]
        String out_dir
    }

    scatter (vcf_info in ptato_vcfs) {
        Pair[String, String] ids = vcf_info.left
        File vcf = vcf_info.right
        
        call PtatoFilter.ptatoIndelFilter {
            input:
                vcf = vcf,
                donor_id = ids.left,
                sample_id = ids.right
        }

        call Htslib.bgzip {
            input:
                file = ptatoIndelFilter.filtered_vcf
        }

        call Htslib.tabix {
            input:
                gz_file = bgzip.gz_file
        }

        # Create output directory structure and copy files
        String output_path = out_dir + "/indels/" + ids.left + "/" + ids.right
        
        call copy_outputs {
            input:
                vcf_gz = bgzip.gz_file,
                vcf_tbi = tabix.tbi_file,
                output_path = output_path
        }
    }

    output {
        Array[Pair[Pair[String, String], Pair[File, File]]] filtered_vcfs = zip(
            ptato_vcfs.left,
            zip(copy_outputs.out_vcf_gz, copy_outputs.out_vcf_tbi)
        )
    }
}

# Task to handle file copying to output directory
task copy_outputs {
    input {
        File vcf_gz
        File vcf_tbi
        String output_path
    }

    command <<<
        mkdir -p ~{output_path}
        cp ~{vcf_gz} ~{output_path}/
        cp ~{vcf_tbi} ~{output_path}/
    >>>

    output {
        File out_vcf_gz = output_path + "/" + basename(vcf_gz)
        File out_vcf_tbi = output_path + "/" + basename(vcf_tbi)
    }

    runtime {
        docker: "ubuntu:latest"
    }
}
