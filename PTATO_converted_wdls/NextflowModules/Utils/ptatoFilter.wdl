version 1.0
workflow PtatoFilterWorkflow {
  input {
    String donor_id
    String sample_id
    File ptaprob_cutoff
    File ptato_vcf
    File ptato_tbi
    File walker_vcf
    File walker_tbi
    File indel_vcf
    File indel_tbi
  }

  call PtatoFilter {
    input:
      donor_id = donor_id,
      sample_id = sample_id,
      ptaprob_cutoff = ptaprob_cutoff,
      ptato_vcf = ptato_vcf,
      ptato_tbi = ptato_tbi,
      walker_vcf = walker_vcf,
      walker_tbi = walker_tbi
  }

  call PtatoIndelFilter {
    input:
      donor_id = donor_id,
      sample_id = sample_id,
      vcf = indel_vcf,
      tbi = indel_tbi
  }

  output {
    File snvs_filtered_vcf = PtatoFilter.ptato_filtered_vcf
    File indels_filtered_vcf = PtatoIndelFilter.ptato_filtered_indel_vcf
  }
}
task PtatoFilter {
  input {
    String donor_id
    String sample_id
    File ptaprob_cutoff
    File ptato_vcf
    File ptato_tbi
    File walker_vcf
    File walker_tbi
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    Rscript ~{default="/scripts/R/ptatoFilter.R"} \
      ~{ptato_vcf} \
      ~{walker_vcf} \
      ~{ptaprob_cutoff} \
      ~{sample_id}.snvs.ptato.filtered.vcf
  >>>

  output {
    File ptato_filtered_vcf = "~{sample_id}.snvs.ptato.filtered.vcf"
  }

  runtime {
    docker: "r-base:latest"
    cpu: 4
    memory: "8 GB"
  }
}
task PtatoIndelFilter {
  input {
    String donor_id
    String sample_id
    File vcf
    File tbi
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    Rscript ~{default="/scripts/R/ptatoIndelFilter.R"} \
      ~{vcf} \
      ~{sample_id}.indels.ptato.filtered.vcf
  >>>

  output {
    File ptato_filtered_indel_vcf = "~{sample_id}.indels.ptato.filtered.vcf"
  }

  runtime {
    docker: "r-base:latest"
    cpu: 4
    memory: "8 GB"
  }
}
