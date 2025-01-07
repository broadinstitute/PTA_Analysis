version 1.0

# Description:
# This workflow generates a genomic context BED file from a VCF file for a given sample.
# The context BED file contains information about genomic regions relevant to further analysis,
# such as indel exclusion or SV detection. The task uses a reference genome to extract context-specific regions.
# The output is a BED file named according to the sample ID.

workflow GetContextWorkflow {
  input {
    String donor_id
    String sample_id
    File vcf
    File tbi
    File ref_genome
  }

  call GetContext {
    input:
      donor_id = donor_id,
      sample_id = sample_id,
      vcf = vcf,
      tbi = tbi,
      ref_genome = ref_genome
  }

  output {
    File bed = GetContext.bed
  }
}

task GetContext {
  input {
    String donor_id
    String sample_id
    File vcf
    File tbi
    File ref_genome
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    R --slave --file=scripts/R/get_context.R --args \
      ~{vcf} \
      ~{sample_id}.context.bed \
      ~{ref_genome}
  >>>

  output {
    File bed = "~{sample_id}.context.bed"
  }

  runtime {
    docker: "r-base:latest"
    cpu: 2
    memory: "4 GB"
  }
}
