version 1.0

# Description:
# This workflow filters Cobalt readcount data for tumor and normal samples.
# The Cobalt readcount file provides coverage data across genomic bins, which is critical for copy number variation (CNV) analysis.
# The task uses provided centromere, cytoband, panel of normals (PoN), and reference genome information.
# The output is a set of filtered readcount files for downstream CNV analysis.

workflow FilterCobaltWorkflow {
  input {
    String donor_id
    String normal_sample_id
    String tumor_sample_id
    File cobalt_tsv_file
    File centromeres
    File cytoband
    File pon  # Panel of normals for filtering
    File ref_genome
    String optional_params = ""  # Optional parameters for filtering
  }

  call FilterCobalt {
    input:
      donor_id = donor_id,
      normal_sample_id = normal_sample_id,
      tumor_sample_id = tumor_sample_id,
      cobalt_tsv_file = cobalt_tsv_file,
      centromeres = centromeres,
      cytoband = cytoband,
      pon = pon,
      ref_genome = ref_genome,
      optional_params = optional_params
  }

  output {
    Array[File] cobalt_filtered_files = FilterCobalt.cobalt_filtered_files
  }
}

task FilterCobalt {
  input {
    String donor_id
    String normal_sample_id
    String tumor_sample_id
    File cobalt_tsv_file
    File centromeres
    File cytoband
    File pon  # Panel of normals file
    File ref_genome
    String optional_params
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    R --slave --file=scripts/R/filter_cobalt.R --args \
      ~{cobalt_tsv_file} \
      ~{centromeres} \
      ~{cytoband} \
      ~{pon} \
      ~{ref_genome} \
      ~{tumor_sample_id} \
      ~{optional_params}
  >>>

  output {
    Array[File] cobalt_filtered_files = glob("~{tumor_sample_id}.readcounts.*")
  }

  runtime {
    docker: "r-base:latest"
    cpu: 4
    memory: "8 GB"
  }
}
