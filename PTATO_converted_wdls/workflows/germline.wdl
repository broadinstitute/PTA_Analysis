version 1.0

import "../tasks/get_gzipped_vcfs.wdl" as gzipped_vcfs
import "../tasks/snpsift.wdl" as snpsift_tasks
import "../tasks/htslib.wdl" as htslib

workflow get_germline_vcfs {
    input {
        # Input files
        Array[File] input_vcfs
        Array[String] bulk_names
        String? germline_vcfs_dir
        String out_dir
    }

    if (defined(germline_vcfs_dir)) {
        call extract_germline_vcf_from_dir {
            input:
                vcfs_dir = select_first([germline_vcfs_dir])
        }

        call gzipped_vcfs.get_gzipped_vcfs {
            input:
                raw_vcfs = extract_germline_vcf_from_dir.vcfs
        }
    }

    if (!defined(germline_vcfs_dir)) {
        scatter (idx in range(length(input_vcfs))) {
            call snpsift_tasks.SnpSift {
                input:
                    vcf = input_vcfs[idx],
                    sample_id = bulk_names[idx]
            }

            call htslib.bgzip {
                input:
                    vcf = SnpSift.filtered_vcf
            }

            call htslib.tabix {
                input:
                    vcf_gz = bgzip.zipped_vcf
            }

            # Copy files to output directory
            call copy_germline_files {
                input:
                    vcf_gz = bgzip.zipped_vcf,
                    vcf_tbi = tabix.vcf_index,
                    donor_id = bulk_names[idx],
                    sample_id = bulk_names[idx],
                    out_dir = out_dir
            }
        }
    }

    output {
        Array[File] germline_vcfs = select_first([
            get_gzipped_vcfs.gzipped_vcfs,
            copy_germline_files.output_vcfs
        ])
        Array[File] germline_vcf_indexes = select_first([
            get_gzipped_vcfs.vcf_indexes,
            copy_germline_files.output_indexes
        ])
    }
}

task extract_germline_vcf_from_dir {
    input {
        String vcfs_dir
    }

    command <<<
        find ~{vcfs_dir} -name "*.vcf" -o -name "*.vcf.gz"
    >>>

    output {
        Array[File] vcfs = read_lines(stdout())
    }

    runtime {
        docker: "ubuntu:latest"
    }
}

task copy_germline_files {
    input {
        File vcf_gz
        File vcf_tbi
        String donor_id
        String sample_id
        String out_dir
    }

    command <<<
        mkdir -p ~{out_dir}/intermediate/germline/~{donor_id}
        cp ~{vcf_gz} ~{out_dir}/intermediate/germline/~{donor_id}/
        cp ~{vcf_tbi} ~{out_dir}/intermediate/germline/~{donor_id}/
    >>>

    output {
        File output_vcfs = "~{out_dir}/intermediate/germline/~{donor_id}/~{basename(vcf_gz)}"
        File output_indexes = "~{out_dir}/intermediate/germline/~{donor_id}/~{basename(vcf_tbi)}"
    }

    runtime {
        docker: "ubuntu:latest"
    }
}