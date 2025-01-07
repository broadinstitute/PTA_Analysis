version 1.0

# Description:
# This workflow filters B-allele frequency (BAF) data from germline VCF files for tumor and normal samples.
# The filtering process uses provided centromere, cytoband, and reference genome information.
# The output is a set of filtered BAF files for downstream structural variant analysis.


workflow FilterBAFWorkflow {
  input {
    String donor_id
    String normal_sample_id
    String tumor_sample_id
    File germline_vcf_file
    File germline_tbi
    File centromeres
    File cytoband
    File ref_genome
    String optional_params = ""  # Optional parameters for filtering
  }

  call FilterBAF {
    input:
      donor_id = donor_id,
      normal_sample_id = normal_sample_id,
      tumor_sample_id = tumor_sample_id,
      germline_vcf_file = germline_vcf_file,
      germline_tbi = germline_tbi,
      centromeres = centromeres,
      cytoband = cytoband,
      ref_genome = ref_genome,
      optional_params = optional_params
  }

  output {
    Array[File] baf_filtered_files = FilterBAF.baf_filtered_files
  }
}

task FilterBAF {
  input {
    String donor_id
    String normal_sample_id
    String tumor_sample_id
    File germline_vcf_file
    File germline_tbi
    File centromeres
    File cytoband
    File ref_genome
    String optional_params
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    R --slave --file=scripts/R/filter_baf.R --args \
      ~{germline_vcf_file} \
      ~{tumor_sample_id} \
      ~{normal_sample_id} \
      ~{centromeres} \
      ~{cytoband} \
      ~{ref_genome} \
      ~{tumor_sample_id} \
      ~{optional_params}
  >>>

  output {
    Array[File] baf_filtered_files = glob("~{tumor_sample_id}.baf.*")
  }

  runtime {
    docker: "r-base:latest"
    cpu: 4
    memory: "8 GB"
  }
}
