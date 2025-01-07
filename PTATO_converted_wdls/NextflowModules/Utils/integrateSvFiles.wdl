version 1.0
workflow IntegrateSvFilesWorkflow {
  input {
    String donor_id
    String normal_sample_id
    String tumor_sample_id
    File baf_filtered_file
    File readcounts_file
    File baf_segments
    File readcounts_segments
    File gripss_filtered_vcf
    String optional_params
  }

  call IntegrateSvFiles {
    input:
      donor_id = donor_id,
      normal_sample_id = normal_sample_id,
      tumor_sample_id = tumor_sample_id,
      baf_filtered_file = baf_filtered_file,
      readcounts_file = readcounts_file,
      baf_segments = baf_segments,
      readcounts_segments = readcounts_segments,
      gripss_filtered_vcf = gripss_filtered_vcf,
      optional_params = optional_params
  }

  output {
    File integrated_sv_files = IntegrateSvFiles.integrated_sv_files
  }
}
task IntegrateSvFiles {
  input {
    String donor_id
    String normal_sample_id
    String tumor_sample_id
    File baf_filtered_file
    File readcounts_file
    File baf_segments
    File readcounts_segments
    File gripss_filtered_vcf
    String optional_params
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    Rscript ~{default="/scripts/R/integrate_sv_files.R"} \
      ~{baf_filtered_file} \
      ~{readcounts_file} \
      ~{baf_segments} \
      ~{readcounts_segments} \
      ~{gripss_filtered_vcf} \
      ~{tumor_sample_id} \
      ~{tumor_sample_id} \
      ~{optional_params}
  >>>

  output {
    File integrated_sv_files = "~{tumor_sample_id}.integrated.*"
  }

  runtime {
    docker: "r-base:latest"
    cpu: 4
    memory: "8 GB"
  }
}
