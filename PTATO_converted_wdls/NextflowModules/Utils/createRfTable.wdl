version 1.0

workflow CreateRfTablesWorkflow {
  input {
    String donor_id
    String sample_id
    File ab_table
    File bed
  }

  call CreateSnvRfTable {
    input:
      donor_id = donor_id,
      sample_id = sample_id,
      ab_table = ab_table,
      bed = bed
  }

  call CreateIndelRfTable {
    input:
      donor_id = donor_id,
      sample_id = sample_id,
      ab_table = ab_table,
      bed = bed
  }

  output {
    File snv_rf_table = CreateSnvRfTable.rf_table
    File indel_rf_table = CreateIndelRfTable.rf_table
  }
}

task CreateSnvRfTable {
  input {
    String donor_id
    String sample_id
    File ab_table
    File bed
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    R --slave --file=scripts/R/create_snv_rf_table.R --args ~{ab_table} ~{bed} ~{donor_id} ~{sample_id}.rftable.rds
  >>>

  output {
    File rf_table = "~{sample_id}.rftable.rds"
  }

  runtime {
    docker: "r-base:latest"
    cpu: 2
    memory: "4 GB"
  }
}

task CreateIndelRfTable {
  input {
    String donor_id
    String sample_id
    File ab_table
    File bed
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    R --slave --file=scripts/R/create_indel_rf_table.R --args ~{ab_table} ~{bed} ~{donor_id} ~{sample_id}.rftable.rds
  >>>

  output {
    File rf_table = "~{sample_id}.rftable.rds"
  }

  runtime {
    docker: "r-base:latest"
    cpu: 2
    memory: "4 GB"
  }
}
