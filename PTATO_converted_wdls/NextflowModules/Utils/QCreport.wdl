version 1.0

workflow QCreportWorkflow {
  input {
    String donor_id
    Array[String] sample_ids
    Array[File] insert_size_metrics_files
    Array[File] wgs_metrics_files
  }

  call QCreport {
    input:
      donor_id = donor_id,
      sample_ids = sample_ids,
      insert_size_metrics_files = insert_size_metrics_files,
      wgs_metrics_files = wgs_metrics_files
  }

  output {
    File qc_report_pdf = QCreport.qc_report_pdf
    File qc_report_txt = QCreport.qc_report_txt
  }
}

task QCreport {
  input {
    String donor_id
    Array[String] sample_ids
    Array[File] insert_size_metrics_files
    Array[File] wgs_metrics_files
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    input_args_1=$(IFS=,; echo "~{sep=',' sample_ids}")
    input_args_2=$(IFS=,; echo "~{sep=',' insert_size_metrics_files}")
    input_args_3=$(IFS=,; echo "~{sep=',' wgs_metrics_files}")

    Rscript ~{default="/scripts/R/PTA_QC_report.R"} \
      ${input_args_1} \
      ${input_args_2} \
      ${input_args_3} \
      ~{donor_id}.qcreport.pdf
  >>>

  output {
    File qc_report_pdf = "~{donor_id}.qcreport.pdf"
    File qc_report_txt = "~{donor_id}.qcreport.txt"
  }

  runtime {
    docker: "r-base:latest"
    cpu: 4
    memory: "8 GB"
  }
}
