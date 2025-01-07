version 1.0
workflow TrainTestRandomForestWorkflow {
  input {
    String label_1
    File rf_table_1
    String label_2
    File rf_table_2
    String train_version
    String donor_id
    String sample_id
    File somatic_vcf
    File somatic_tbi
    File rf_model_rds_snv
    File rf_model_rds_indel
  }

  call TrainSnvRf {
    input:
      label_1 = label_1,
      rf_table_1 = rf_table_1,
      label_2 = label_2,
      rf_table_2 = rf_table_2,
      train_version = train_version
  }

  call TrainIndelRf {
    input:
      label_1 = label_1,
      rf_table_1 = rf_table_1,
      label_2 = label_2,
      rf_table_2 = rf_table_2,
      train_version = train_version
  }

  call TestSnvRf {
    input:
      donor_id = donor_id,
      sample_id = sample_id,
      somatic_vcf = somatic_vcf,
      somatic_tbi = somatic_tbi,
      rf_table = rf_table_1,
      rf_model_rds = rf_model_rds_snv
  }

  call TestIndelRf {
    input:
      donor_id = donor_id,
      sample_id = sample_id,
      somatic_vcf = somatic_vcf,
      somatic_tbi = somatic_tbi,
      rf_table = rf_table_2,
      rf_model_rds = rf_model_rds_indel
  }


task TrainSnvRf {
  input {
    String label_1
    File rf_table_1
    String label_2
    File rf_table_2
    String train_version
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    Rscript ~{default="/scripts/R/train_randomforest.R"} \
      ~{label_1} ~{rf_table_1} ~{label_2} ~{rf_table_2} ~{train_version}
  >>>

  output {
    File confusion_txt = glob("randomforest*confusion.txt")[0]
    File importance_txt = glob("randomforest*importance.txt")[0]
    File rdata_file = glob("randomforest*.Rdata")[0]
    File rds_file = glob("randomforest*.rds")[0]
  }

  runtime {
    docker: "r-base:latest"
    cpu: 4
    memory: "8 GB"
  }
}

task TrainIndelRf {
  input {
    String label_1
    File rf_table_1
    String label_2
    File rf_table_2
    String train_version
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    Rscript ~{default="/scripts/R/train_randomforest.R"} \
      ~{label_1} ~{rf_table_1} ~{label_2} ~{rf_table_2} ~{train_version}
  >>>

  output {
    File confusion_txt = glob("randomforest*confusion.txt")[0]
    File importance_txt = glob("randomforest*importance.txt")[0]
    File rdata_file = glob("randomforest*.Rdata")[0]
    File rds_file = glob("randomforest*.rds")[0]
  }

  runtime {
    docker: "r-base:latest"
    cpu: 4
    memory: "8 GB"
  }
}
task TestSnvRf {
  input {
    String donor_id
    String sample_id
    File somatic_vcf
    File somatic_tbi
    File rf_table
    File rf_model_rds
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    Rscript ~{default="/scripts/R/test_randomforest.R"} \
      ~{rf_model_rds} ~{rf_table} ~{somatic_vcf} ~{sample_id}.snvs.ptato.vcf
  >>>

  output {
    File ptato_vcf = "~{sample_id}.snvs.ptato.vcf"
  }

  runtime {
    docker: "r-base:latest"
    cpu: 4
    memory: "8 GB"
  }
}
task TestIndelRf {
  input {
    String donor_id
    String sample_id
    File somatic_vcf
    File somatic_tbi
    File rf_table
    File rf_model_rds
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    Rscript ~{default="/scripts/R/test_randomforest.R"} \
      ~{rf_model_rds} ~{rf_table} ~{somatic_vcf} ~{sample_id}.indels.ptato.vcf
  >>>

  output {
    File ptato_vcf = "~{sample_id}.indels.ptato.vcf"
  }

  runtime {
    docker: "r-base:latest"
    cpu: 4
    memory: "8 GB"
  }
}
