version 1.0
workflow PtatoCutoffWorkflow {
  input {
    String donor_id
    String sample_id
    File ptato_vcf
    File ptato_tbi
    File walker_vcf
    File walker_tbi
    String ref_genome
    String optional_params
  }

  call PtatoCutoff {
    input:
      donor_id = donor_id,
      sample_id = sample_id,
      ptato_vcf = ptato_vcf,
      ptato_tbi = ptato_tbi,
      walker_vcf = walker_vcf,
      walker_tbi = walker_tbi,
      ref_genome = ref_genome,
      optional_params = optional_params
  }

  output {
    File ptato_table = PtatoCutoff.ptato_table
    File ptaprob_cutoff = PtatoCutoff.ptaprob_cutoff
  }
}
task PtatoCutoff {
  input {
    String donor_id
    String sample_id
    File ptato_vcf
    File ptato_tbi
    File walker_vcf
    File walker_tbi
    String ref_genome
    String optional_params
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    Rscript ~{default="/scripts/R/ptatoCutoff.R"} \
      ~{ptato_vcf} \
      ~{walker_vcf} \
      ~{ref_genome} \
      ~{optional_params} \
      ~{sample_id}.ptatotable.txt > ~{sample_id}.ptaprobcutoff.txt
  >>>

  output {
    File ptato_table = "~{sample_id}.ptatotable.txt"
    File ptaprob_cutoff = "~{sample_id}.ptaprobcutoff.txt"
  }

  runtime {
    docker: "r-base:latest"
    cpu: 4
    memory: "8 GB"
  }
}
