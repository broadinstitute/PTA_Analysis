version 1.0

workflow ExcludeIndelsWorkflow {
  input {
    String donor_id
    String sample_id
    File indels_vcf
    File indels_tbi
    File context_bed
    File exclude_vcf
  }

  call ExcludeIndels {
    input:
      donor_id = donor_id,
      sample_id = sample_id,
      indels_vcf = indels_vcf,
      indels_tbi = indels_tbi,
      context_bed = context_bed,
      exclude_vcf = exclude_vcf
  }

  output {
    File ptato_vcf = ExcludeIndels.ptato_vcf
  }
}

task ExcludeIndels {
  input {
    String donor_id
    String sample_id
    File indels_vcf
    File indels_tbi
    File context_bed
    File exclude_vcf
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    R --slave --file=scripts/R/excludeIndels.R --args \
      ~{indels_vcf} \
      ~{context_bed} \
      ~{exclude_vcf} \
      ~{sample_id}.indels.ptato.vcf
  >>>

  output {
    File ptato_vcf = "~{sample_id}.indels.ptato.vcf"
  }

  runtime {
    docker: "r-base:latest"
    cpu: 2
    memory: "4 GB"
  }
}
