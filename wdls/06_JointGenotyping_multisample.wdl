version 1.0

# The Broad Institute Software License
# Copyright (c) 2025 The Broad Institute, Inc.
# All rights reserved.
#
# This software is made available under the Broad Institute Software License.
# You may use, copy, modify, and distribute this software under the terms of the license.
# 
# For details, please see https://software.broadinstitute.org/gatk/

# Workflow: GenotypeGVCFs
# Description: Performs joint genotyping on a single multi-sample GVCF file (e.g., from CombineGVCFs)
#              using GATK's GenotypeGVCFs following Broad best practices.
# Author: Shadi Zaheri
# Date: 2025-05-10

workflow GenotypeGVCFs {
  parameter_meta {
    reference_fasta: "Path to the reference genome in FASTA format."
    reference_fasta_index: "Index for the reference FASTA (.fai)."
    reference_dict: "Sequence dictionary (.dict) for the reference."
    combined_gvcf: "Input GVCF file from CombineGVCFs."
    combined_gvcf_index: "Index (.tbi) for the combined GVCF."
    output_vcf_basename: "Base name for the final genotyped VCF."
    dbsnp_vcf: "Optional dbSNP VCF file for annotation."
    dbsnp_vcf_index: "Index file (.tbi) for dbSNP VCF."
    interval_list: "Optional BED or interval list file for targeted genotyping."
    preemptible_tries: "Number of preemptible retries allowed."
    memory_gb: "Memory allocated to the task (in GB)."
    disk_gb: "Disk size allocated (in GB)."
    cpu: "Number of CPU cores used."
    gatk_docker: "Docker image URI for GATK (default: broadinstitute/gatk:4.6.1.0)."
  }

  input {
    File reference_fasta
    File reference_fasta_index
    File reference_dict
    File combined_gvcf
    File combined_gvcf_index
    String output_vcf_basename
    File? dbsnp_vcf
    File? dbsnp_vcf_index
    File? interval_list
    Int preemptible_tries = 1
    Int memory_gb = 32
    Int disk_gb = 200
    Int cpu = 4
    String gatk_docker = "broadinstitute/gatk:4.6.1.0"
  }

  call GenotypeGVCFsTask {
    input:
      reference_fasta = reference_fasta,
      reference_fasta_index = reference_fasta_index,
      reference_dict = reference_dict,
      combined_gvcf = combined_gvcf,
      combined_gvcf_index = combined_gvcf_index,
      output_vcf_basename = output_vcf_basename,
      dbsnp_vcf = dbsnp_vcf,
      dbsnp_vcf_index = dbsnp_vcf_index,
      interval_list = interval_list,
      preemptible_tries = preemptible_tries,
      memory_gb = memory_gb,
      disk_gb = disk_gb,
      cpu = cpu,
      gatk_docker = gatk_docker
  }

  output {
    File genotyped_vcf = GenotypeGVCFsTask.genotyped_vcf
    File genotyped_vcf_index = GenotypeGVCFsTask.genotyped_vcf_index
  }
}

task GenotypeGVCFsTask {
  parameter_meta {
    reference_fasta: "Reference genome FASTA."
    reference_fasta_index: "Index file for reference FASTA."
    reference_dict: "Sequence dictionary for reference FASTA."
    combined_gvcf: "Combined GVCF input."
    combined_gvcf_index: "Index (.tbi) for combined GVCF."
    output_vcf_basename: "Base name for output VCF."
    dbsnp_vcf: "Optional dbSNP file for annotation."
    dbsnp_vcf_index: "Index for dbSNP file."
    interval_list: "Optional list of intervals (e.g., BED)."
    preemptible_tries: "Number of retries on preemptible VMs."
    memory_gb: "Memory allocated (in GB)."
    disk_gb: "Disk size allocated (in GB)."
    cpu: "Number of CPUs."
    gatk_docker: "Docker image to use (GATK 4.x expected)."
  }

  input {
    File reference_fasta
    File reference_fasta_index
    File reference_dict
    File combined_gvcf
    File combined_gvcf_index
    String output_vcf_basename
    File? dbsnp_vcf
    File? dbsnp_vcf_index
    File? interval_list
    Int preemptible_tries
    Int memory_gb
    Int disk_gb
    Int cpu
    String gatk_docker
  }

  command {
    gatk --java-options "-Xmx~{memory_gb}G" GenotypeGVCFs \
      -R ~{reference_fasta} \
      --variant ~{combined_gvcf} \
      -O ~{output_vcf_basename}.vcf.gz \
      ~{if defined(dbsnp_vcf) then "--dbsnp " + dbsnp_vcf else ""} \
      ~{if defined(interval_list) then "-L " + interval_list else ""}
  }

  runtime {
    docker: gatk_docker
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: preemptible_tries
    cpu: cpu
  }

  output {
    File genotyped_vcf = "~{output_vcf_basename}.vcf.gz"
    File genotyped_vcf_index = "~{output_vcf_basename}.vcf.gz.tbi"
  }
}
