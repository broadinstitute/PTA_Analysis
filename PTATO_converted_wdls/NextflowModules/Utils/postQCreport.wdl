version 1.0
workflow PostQCreportWorkflow {
  input {
    String donor_id
    Array[String] sample_ids
    File ptato_vcf
    File ptato_tbi
    File ptato_filt_vcf
    File ptato_filt_tbi
    File ptato_table
    File walker_vcf
    File walker_tbi
    File autosomal_callable_files
  }

  call PostQCreport {
    input:
      donor_id = donor_id,
      sample_ids = sample_ids,
      ptato_vcf = ptato_vcf,
      ptato_tbi = ptato_tbi,
      ptato_filt_vcf = ptato_filt_vcf,
      ptato_filt_tbi = ptato_filt_tbi,
      ptato_table = ptato_table,
      walker_vcf = walker_vcf,
      walker_tbi = walker_tbi,
      autosomal_callable_files = autosomal_callable_files
  }

  output {
    File postqc_report_pdf = PostQCreport.postqc_report_pdf
    File postqc_report_txt = PostQCreport.postqc_report_txt
  }
}

task PostQCreport {
  input {
    String donor_id
    Array[String] sample_ids
    File ptato_vcf
    File ptato_tbi
    File ptato_filt_vcf
    File ptato_filt_tbi
    File ptato_table
    File walker_vcf
    File walker_tbi
    File autosomal_callable_files
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    input_args_1=$(IFS=,; echo "~{sep=',' sample_ids}")
    input_args_2="~{ptato_vcf}"
    input_args_3="~{ptato_filt_vcf}"
    input_args_4="~{walker_vcf}"
    input_args_5="~{ptato_table}"
    input_args_6="~{autosomal_callable_files}"

    Rscript ~{default="/scripts/R/PostPTATO_QC_report.R"} \
      ${input_args_1} \
      ${input_args_2} \
      ${input_args_3} \
      ${input_args_4} \
      ${input_args_5} \
      ${input_args_6} \
      ~{donor_id}.postqcreport
  >>>

  output {
    File postqc_report_pdf = "~{donor_id}.postqcreport.pdf"
    File postqc_report_txt = "~{donor_id}.postqcreport.txt"
  }

  runtime {
    docker: "r-base:latest"
    cpu: 4
    memory: "8GB"
  }
}
