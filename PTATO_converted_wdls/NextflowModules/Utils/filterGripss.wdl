version 1.0

# Description:
# This workflow filters structural variants (SVs) detected by GRIPSS for tumor and normal samples.
# GRIPSS provides somatic SV calls, which are further refined using B-allele frequency (BAF) data,
# Cobalt readcount data, and a panel of normals (PoN) to remove false positives.
# The output is a set of filtered SV files for downstream analysis.

workflow FilterGripssWorkflow {
  input {
    String donor_id
    String normal_sample_id
    String tumor_sample_id
    File gripss_somatic_filtered_vcf
    File gripss_somatic_filtered_tbi
    File baf_filtered_file
    File cobalt_filtered_readcounts_file
    File pon  # Panel of normals for filtering
    String optional_params = ""  # Optional parameters for filtering
  }

  call FilterGripss {
    input:
      donor_id = donor_id,
      normal_sample_id = normal_sample_id,
      tumor_sample_id = tumor_sample_id,
      gripss_somatic_filtered_vcf = gripss_somatic_filtered_vcf,
      gripss_somatic_filtered_tbi = gripss_somatic_filtered_tbi,
      baf_filtered_file = baf_filtered_file,
      cobalt_filtered_readcounts_file = cobalt_filtered_readcounts_file,
      pon = pon,
      optional_params = optional_params
  }

  output {
    Array[File] gripss_filtered_files = FilterGripss.gripss_filtered_files
  }
}

task FilterGripss {
  input {
    String donor_id
    String normal_sample_id
    String tumor_sample_id
    File gripss_somatic_filtered_vcf
    File gripss_somatic_filtered_tbi
    File baf_filtered_file
    File cobalt_filtered_readcounts_file
    File pon  # Panel of normals file
    String optional_params
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    R --slave --file=scripts/R/filter_gripss.R --args \
      ~{gripss_somatic_filtered_vcf} \
      ~{baf_filtered_file} \
      ~{cobalt_filtered_readcounts_file} \
      ~{pon} \
      ~{tumor_sample_id} \
      ~{optional_params}
  >>>

  output {
    Array[File] gripss_filtered_files = glob("~{tumor_sample_id}.svs.*")
  }

  runtime {
    docker: "r-base:latest"
    cpu: 4
    memory: "8 GB"
  }
}
