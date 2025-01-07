version 1.0

workflow CreateSvPlotsWorkflow {
  input {
    String donor_id
    String normal_sample_id
    String tumor_sample_id
    File readcounts_100kb_file
    File readcounts_1mb_file
    File readcounts_segments_file
    File baf_binned_100kb_file
    File baf_segments_file
    File cnv_file
    String optional_params = ""  # Optional parameters for the script
  }

  call CreateSvPlots {
    input:
      donor_id = donor_id,
      normal_sample_id = normal_sample_id,
      tumor_sample_id = tumor_sample_id,
      readcounts_100kb_file = readcounts_100kb_file,
      readcounts_1mb_file = readcounts_1mb_file,
      readcounts_segments_file = readcounts_segments_file,
      baf_binned_100kb_file = baf_binned_100kb_file,
      baf_segments_file = baf_segments_file,
      cnv_file = cnv_file,
      optional_params = optional_params
  }

  output {
    Array[File] sv_plots = CreateSvPlots.sv_plots
  }
}

task CreateSvPlots {
  input {
    String donor_id
    String normal_sample_id
    String tumor_sample_id
    File readcounts_100kb_file
    File readcounts_1mb_file
    File readcounts_segments_file
    File baf_binned_100kb_file
    File baf_segments_file
    File cnv_file
    String optional_params
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    R --slave --file=scripts/R/create_sv_plots.R --args \
      ~{readcounts_100kb_file} \
      ~{readcounts_1mb_file} \
      ~{readcounts_segments_file} \
      ~{baf_binned_100kb_file} \
      ~{baf_segments_file} \
      ~{cnv_file} \
      ~{tumor_sample_id} \
      ~{tumor_sample_id} \
      ~{optional_params}
  >>>

  output {
    Array[File] sv_plots = glob("~{tumor_sample_id}.*.p*")
  }

  runtime {
    docker: "r-base:latest"
    cpu: 4
    memory: "8 GB"
  }
}
