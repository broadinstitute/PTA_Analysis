version 1.0

# Workflow: BaseQualityScoreRecalibration
# Author: Shadi Zaheri
# Date: 2025-01-26
# Description: Performs Base Quality Score Recalibration (BQSR) on WGS data using GATK BaseRecalibrator and ApplyBQSR. 
# The input BAM should be the output from MarkDuplicates, and the workflow assumes known variant databases are available.

workflow BaseQualityScoreRecalibration {
  parameter_meta {
    input_bam: "Path to the input BAM file. Must be the output BAM from MarkDuplicates."
    input_bam_index: "Path to the index file for the input BAM."
    reference_fasta: "Path to the reference genome in FASTA format."
    reference_fasta_index: "Path to the index file for the reference FASTA."
    reference_dict: "Path to the sequence dictionary for the reference FASTA."
    known_sites_vcf: "Array of known variant sites for recalibration (e.g., dbSNP, Mills, 1000G)."
    known_sites_index: "Array of index files (.tbi) for the known variant sites VCFs."
    recalibration_report_filename: "Filename for the recalibration report (default includes sample name)."
    sample_name: "Unique sample name used to name output files."
    output_bam_basename: "Base name for the recalibrated BAM output (default includes sample name)."
    preemptible_tries: "Number of preemptible retries allowed for each task."
    memory_gb: "Memory allocated for each task in gigabytes."
    disk_gb: "Disk space allocated for each task in gigabytes."
    cpu: "Number of CPU cores allocated for each task."
    gatk_docker: "Docker image for GATK tools (default: broadinstitute/gatk:4.6.1.0)."
  }

  input {
    File input_bam
    File input_bam_index
    File reference_fasta
    File reference_fasta_index
    File reference_dict
    Array[File] known_sites_vcf
    Array[File] known_sites_index
    String sample_name
    String recalibration_report_filename = "~{sample_name}_recal_data.table"
    String output_bam_basename = "~{sample_name}_recalibrated"
    Int preemptible_tries = 1
    Int memory_gb = 16
    Int disk_gb = 200
    Int cpu = 4
    String gatk_docker = "broadinstitute/gatk:4.6.1.0"
  }

  call BaseRecalibrator {
    input:
      input_bam = input_bam,
      input_bam_index = input_bam_index,
      reference_fasta = reference_fasta,
      reference_fasta_index = reference_fasta_index,
      reference_dict = reference_dict,
      known_sites_vcf = known_sites_vcf,
      known_sites_index = known_sites_index,
      recalibration_report_filename = recalibration_report_filename,
      preemptible_tries = preemptible_tries,
      memory_gb = memory_gb,
      disk_gb = disk_gb,
      cpu = cpu,
      gatk_docker = gatk_docker
  }

  call ApplyBQSR {
    input:
      input_bam = input_bam,
      input_bam_index = input_bam_index,
      reference_fasta = reference_fasta,
      reference_fasta_index = reference_fasta_index,
      reference_dict = reference_dict,
      recalibration_report = BaseRecalibrator.recalibration_report,
      output_bam_basename = output_bam_basename,
      preemptible_tries = preemptible_tries,
      memory_gb = memory_gb,
      disk_gb = disk_gb,
      cpu = cpu,
      gatk_docker = gatk_docker
  }

  output {
    File recalibrated_bam = ApplyBQSR.recalibrated_bam
    File recalibrated_bam_index = ApplyBQSR.recalibrated_bam_index
    File recalibration_report = BaseRecalibrator.recalibration_report
  }
}


task BaseRecalibrator {
  input {
    File input_bam
    File input_bam_index
    File reference_fasta
    File reference_fasta_index
    File reference_dict
    Array[File] known_sites_vcf
    Array[File] known_sites_index
    String recalibration_report_filename
    Int preemptible_tries
    Int memory_gb
    Int disk_gb
    Int cpu
    String gatk_docker
  }

  command {
    set -euxo pipefail

    # Flatten known-sites into paired arguments
    gatk --java-options "-Xmx~{memory_gb}G" BaseRecalibrator \
      -R ~{reference_fasta} \
      -I ~{input_bam} \
      ~{sep=" " prefix("--known-sites ", known_sites_vcf)} \
      -O ~{recalibration_report_filename}
  }

  runtime {
    docker: gatk_docker
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_gb} HDD"
    cpu: cpu
    preemptible: preemptible_tries
  }

  output {
    File recalibration_report = recalibration_report_filename
  }

  meta {
    description: "Runs GATK BaseRecalibrator on the input BAM file, using known sites for recalibration."
  }
}

# Task: ApplyBQSR
task ApplyBQSR {
  input {
    File input_bam
    File input_bam_index
    File reference_fasta
    File reference_fasta_index
    File reference_dict
    File recalibration_report
    String output_bam_basename
    Int preemptible_tries
    Int memory_gb
    Int disk_gb
    Int cpu
    String gatk_docker
  }

  command {
    set -euxo pipefail

    gatk --java-options "-Xmx~{memory_gb}G" ApplyBQSR \
      -R ~{reference_fasta} \
      --reference-index ~{reference_fasta_index} \
      --reference-dict ~{reference_dict} \
      -I ~{input_bam} \
      --bqsr-recal-file ~{recalibration_report} \
      -O ~{output_bam_basename}.bam
  }

  runtime {
    docker: gatk_docker
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_gb} HDD"
    cpu: cpu
    preemptible: preemptible_tries
  }

  output {
    File recalibrated_bam = "~{output_bam_basename}.bam"
    File recalibrated_bam_index = "~{output_bam_basename}.bam.bai"
  }
}
