version 1.0

# Copyright (c) 2025 Broad Institute, Inc.
# All rights reserved.
#
# This software is made available under the Broad Institute Software License Agreement.
# For full details, see: https://software.broadinstitute.org/software/license

workflow CombineGVCFs {
  meta {
    author: "Shadi Zaheri"
    description: "Combines multiple GVCF files into a single multi-sample GVCF for joint genotyping using GATK CombineGVCFs."
  }

  parameter_meta {
    reference_fasta: "Path to the reference genome in FASTA format."
    reference_fasta_index: "Path to the index file for the reference FASTA."
    reference_dict: "Path to the sequence dictionary for the reference FASTA."
    input_gvcfs: "Array of input GVCF files to be merged."
    input_gvcf_indexes: "Array of index files (.tbi) for the input GVCFs."
    output_gvcf_basename: "Base name for the combined output GVCF file."
    interval_file: "Optional: Path to an interval list file (e.g., BED or Picard format) for targeted merging."
    preemptible_tries: "Number of preemptible retries allowed for each task."
    memory_gb: "Memory allocated for each task in gigabytes."
    disk_gb: "Disk space allocated for each task in gigabytes."
    cpu: "Number of CPU cores allocated for each task."
    gatk_docker: "Docker image for GATK tools (default: broadinstitute/gatk:4.6.1.0)."
  }

  input {
    File reference_fasta
    File reference_fasta_index
    File reference_dict
    Array[File] input_gvcfs
    Array[File] input_gvcf_indexes
    String output_gvcf_basename
    File? interval_file  # Optional
    Int preemptible_tries = 1
    Int memory_gb = 32
    Int disk_gb = 200
    Int cpu = 4
    String gatk_docker = "broadinstitute/gatk:4.6.1.0"
  }

  call CombineGVCFsTask {
    input:
      reference_fasta = reference_fasta,
      reference_fasta_index = reference_fasta_index,
      reference_dict = reference_dict,
      input_gvcfs = input_gvcfs,
      input_gvcf_indexes = input_gvcf_indexes,
      output_gvcf_basename = output_gvcf_basename,
      interval_file = interval_file,
      preemptible_tries = preemptible_tries,
      memory_gb = memory_gb,
      disk_gb = disk_gb,
      cpu = cpu,
      gatk_docker = gatk_docker
  }

  output {
    File combined_gvcf = CombineGVCFsTask.combined_gvcf
    File combined_gvcf_index = CombineGVCFsTask.combined_gvcf_index
  }
}

# Task: CombineGVCFs
task CombineGVCFsTask {
  meta {
    author: "Shadi Zaheri"
    description: "Uses GATK CombineGVCFs to merge multiple per-sample GVCFs into a single GVCF file."
  }

  parameter_meta {
    reference_fasta: "Path to the reference genome in FASTA format."
    reference_fasta_index: "Path to the index file for the reference FASTA."
    reference_dict: "Path to the sequence dictionary for the reference FASTA."
    input_gvcfs: "Array of input GVCF files to be merged."
    input_gvcf_indexes: "Array of index files (.tbi) for the input GVCFs."
    output_gvcf_basename: "Base name for the combined output GVCF file."
    interval_file: "Optional: Path to an interval list file (e.g., BED or Picard format) for targeted merging."
    preemptible_tries: "Number of preemptible retries allowed for each task."
    memory_gb: "Memory allocated for each task in gigabytes."
    disk_gb: "Disk space allocated for each task in gigabytes."
    cpu: "Number of CPU cores allocated for each task."
    gatk_docker: "Docker image for GATK tools (default: broadinstitute/gatk:4.6.1.0)."
  }

  input {
    File reference_fasta
    File reference_fasta_index
    File reference_dict
    Array[File] input_gvcfs
    Array[File] input_gvcf_indexes
    String output_gvcf_basename
    File? interval_file  # Optional
    Int preemptible_tries
    Int memory_gb
    Int disk_gb
    Int cpu
    String gatk_docker
  }

  command <<<
    gatk --java-options "-Xmx~{memory_gb}G" CombineGVCFs \
      -R ~{reference_fasta} \
      ~{sep=' ' input_gvcfs} \
      -O ~{output_gvcf_basename}.g.vcf.gz \
      ~{if defined(interval_file) then "-L " + interval_file else ""}
  >>>

  runtime {
    docker: gatk_docker
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: preemptible_tries
    cpu: cpu
  }

  output {
    File combined_gvcf = "~{output_gvcf_basename}.g.vcf.gz"
    File combined_gvcf_index = "~{output_gvcf_basename}.g.vcf.gz.tbi"
  }
}
