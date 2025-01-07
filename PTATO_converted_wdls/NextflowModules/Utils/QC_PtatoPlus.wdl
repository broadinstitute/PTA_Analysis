version 1.0
workflow QC_PtatoPlusWorkflow {
  input {
    String donor_id
    Array[String] sample_ids
    Array[File] PTATO_files
    Array[File] PTATO_filt_files
  }

  call QC_PtatoPlus {
    input:
      donor_id = donor_id,
      sample_ids = sample_ids,
      PTATO_files = PTATO_files,
      PTATO_filt_files = PTATO_filt_files
  }

  output {
    File postPTATO_QC_pdf = QC_PtatoPlus.postPTATO_QC_pdf
    File postPTATO_QC_txt = QC_PtatoPlus.postPTATO_QC_txt
  }
}
task QC_PtatoPlus {
  input {
    String donor_id
    Array[String] sample_ids
    Array[File] PTATO_files
    Array[File] PTATO_filt_files
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    input_args_1=$(IFS=,; echo "~{sep=',' sample_ids}")
    input_args_2=$(IFS=,; echo "~{sep=',' PTATO_files}")
    input_args_3=$(IFS=,; echo "~{sep=',' PTATO_filt_files}")

    Rscript ~{default="/scripts/R/PTATO_Report.R"} \
      ${input_args_1} \
      ${input_args_2} \
      ${input_args_3} \
      ~{donor_id}.postPTATO_QC
  >>>

  output {
    File postPTATO_QC_pdf = "~{donor_id}.postPTATO_QC.pdf"
    File postPTATO_QC_txt = "~{donor_id}.postPTATO_QC.txt"
  }

  runtime {
    docker: "r-base:latest"
    cpu: 4
    memory: "8 GB"
  }
}
