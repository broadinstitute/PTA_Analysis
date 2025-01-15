version 1.0

workflow BaseQualityScoreRecalibration {
  input {
    File input_bam
    File input_bam_index
    File reference_fasta
    File reference_fasta_index
    File reference_dict
    Array[File] known_sites
    String recalibration_report_filename
    String output_bam_basename
    Int preemptible_tries = 1
    Int memory_gb = 8
    Int disk_gb = 50
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.3.0.0"
  }

  call BaseRecalibrator {
    input:
      input_bam = input_bam,
      input_bam_index = input_bam_index,
      reference_fasta = reference_fasta,
      reference_fasta_index = reference_fasta_index,
      reference_dict = reference_dict,
      known_sites = known_sites,
      recalibration_report_filename = recalibration_report_filename,
      preemptible_tries = preemptible_tries,
      memory_gb = memory_gb,
      disk_gb = disk_gb,
      gatk_docker = gatk_docker
  }

  call ApplyBQSR {
    input:
      input_bam = input_bam,
      input_bam_index = input_bam_index,
      reference_fasta = reference_fasta,
      reference_fasta_index = reference_fasta_index,
      recalibration_report = BaseRecalibrator.recalibration_report,
      output_bam_basename = output_bam_basename,
      preemptible_tries = preemptible_tries,
      memory_gb = memory_gb,
      disk_gb = disk_gb,
      gatk_docker = gatk_docker
  }

  output {
    File recalibrated_bam = ApplyBQSR.recalibrated_bam
    File recalibrated_bam_index = ApplyBQSR.recalibrated_bam_index
  }
}

### Task for BaseRecalibrator
task BaseRecalibrator {
  input {
    File input_bam
    File input_bam_index
    File reference_fasta
    File reference_fasta_index
    File reference_dict
    Array[File] known_sites
    String recalibration_report_filename
    Int preemptible_tries
    Int memory_gb
    Int disk_gb
    String gatk_docker
  }

  command {
    gatk --java-options "-Xmx~{memory_gb}G" BaseRecalibrator \
      -R ~{reference_fasta} \
      -I ~{input_bam} \
      --known-sites ~{sep=" --known-sites " known_sites} \
      -O ~{recalibration_report_filename}
  }

  runtime {
    docker: gatk_docker
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: preemptible_tries
  }

  output {
    File recalibration_report = recalibration_report_filename
  }
}

### Task for ApplyBQSR
task ApplyBQSR {
  input {
    File input_bam
    File input_bam_index
    File reference_fasta
    File reference_fasta_index
    File recalibration_report
    String output_bam_basename
    Int preemptible_tries
    Int memory_gb
    Int disk_gb
    String gatk_docker
  }

  command {
    gatk --java-options "-Xmx~{memory_gb}G" ApplyBQSR \
      -R ~{reference_fasta} \
      -I ~{input_bam} \
      --bqsr-recal-file ~{recalibration_report} \
      -O ~{output_bam_basename}.bam
  }

  runtime {
    docker: gatk_docker
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: preemptible_tries
  }

  output {
    File recalibrated_bam = "~{output_bam_basename}.bam"
    File recalibrated_bam_index = "~{output_bam_basename}.bam.bai"
  }
}
