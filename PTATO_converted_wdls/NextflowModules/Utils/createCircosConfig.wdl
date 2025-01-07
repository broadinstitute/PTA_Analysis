version 1.0

workflow CreateCircosConfigWorkflow {
  input {
    String donor_id
    String normal_sample_id
    String tumor_sample_id
    File cnv_file
    File readcounts_1mb_file
    File baf_1mb_file
    File readcounts_segments_txt_file
    File baf_segments_txt_file
    File sv_vcf
    String optional_params = ""  # Optional parameters for the script
  }

  call CreateCircosConfig {
    input:
      donor_id = donor_id,
      normal_sample_id = normal_sample_id,
      tumor_sample_id = tumor_sample_id,
      cnv_file = cnv_file,
      readcounts_1mb_file = readcounts_1mb_file,
      baf_1mb_file = baf_1mb_file,
      readcounts_segments_txt_file = readcounts_segments_txt_file,
      baf_segments_txt_file = baf_segments_txt_file,
      sv_vcf = sv_vcf,
      optional_params = optional_params
  }

  output {
    Array[File] circos_config_files = CreateCircosConfig.circos_config_files
  }
}

task CreateCircosConfig {
  input {
    String donor_id
    String normal_sample_id
    String tumor_sample_id
    File cnv_file
    File readcounts_1mb_file
    File baf_1mb_file
    File readcounts_segments_txt_file
    File baf_segments_txt_file
    File sv_vcf
    String optional_params
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    R --slave --file=scripts/R/create_circos_config.R --args \
      ~{cnv_file} \
      ~{readcounts_1mb_file} \
      ~{baf_1mb_file} \
      ~{readcounts_segments_txt_file} \
      ~{baf_segments_txt_file} \
      ~{sv_vcf} \
      ~{tumor_sample_id} \
      ~{optional_params}
  >>>

  output {
    Array[File] circos_config_files = glob("~{tumor_sample_id}.circos*")
  }

  runtime {
    docker: "r-base:latest"
    cpu: 2
    memory: "8 GB"
  }
}
