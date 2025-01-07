version 1.0

workflow MergePtatoVcfsWorkflow {
  input {
    String donor_id
    String sample_id
    File input_vcf
    File input_tbi
    Array[String] ptato_snvs_sample_ids
    Array[File] ptato_snvs_vcfs
    Array[File] ptato_snvs_tbis
    Array[String] ptato_indels_sample_ids
    Array[File] ptato_indels_vcfs
    Array[File] ptato_indels_tbis
  }

  call MergePtatoVcfs {
    input:
      donor_id = donor_id,
      sample_id = sample_id,
      input_vcf = input_vcf,
      input_tbi = input_tbi,
      ptato_snvs_sample_ids = ptato_snvs_sample_ids,
      ptato_snvs_vcfs = ptato_snvs_vcfs,
      ptato_snvs_tbis = ptato_snvs_tbis,
      ptato_indels_sample_ids = ptato_indels_sample_ids,
      ptato_indels_vcfs = ptato_indels_vcfs,
      ptato_indels_tbis = ptato_indels_tbis
  }

  output {
    File ptato_intersect_vcf = MergePtatoVcfs.ptato_intersect_vcf
  }
}


task MergePtatoVcfs {
  input {
    String donor_id
    String sample_id
    File input_vcf
    File input_tbi
    Array[String] ptato_snvs_sample_ids
    Array[File] ptato_snvs_vcfs
    Array[File] ptato_snvs_tbis
    Array[String] ptato_indels_sample_ids
    Array[File] ptato_indels_vcfs
    Array[File] ptato_indels_tbis
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    ptato_vcfs=$(IFS=,; echo "~{sep=',' ptato_snvs_vcfs},~{sep=',' ptato_indels_vcfs}")

    Rscript ~{default="/scripts/R/merge_ptato_vcfs.R"} \
      ~{input_vcf} \
      ${ptato_vcfs} \
      ~{donor_id}.ptato.merged.vcf
  >>>

  output {
    File ptato_intersect_vcf = "~{donor_id}.ptato.merged.vcf"
  }

  runtime {
    docker: "r-base:latest"
    cpu: 4
    memory: "8 GB"
  }
}
