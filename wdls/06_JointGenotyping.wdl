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
# Date: 2025-02-24
# Description: This workflow performs joint genotyping on multiple GVCFs using GATK's GenotypeGVCFs.
#              It follows the settings from the reference paper and Broad's best practices.
#              The workflow allows optional use of a known dbSNP file and an interval list.


workflow GenotypeGVCFs {
  input {
    File reference_fasta
    File reference_fasta_index
    File reference_dict
    Array[File] input_gvcfs
    Array[File] input_gvcf_indexes
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
      input_gvcfs = input_gvcfs,
      input_gvcf_indexes = input_gvcf_indexes,
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
  input {
    File reference_fasta
    File reference_fasta_index
    File reference_dict
    Array[File] input_gvcfs
    Array[File] input_gvcf_indexes
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
      ~{sep=' ' prefix("--variant ", input_gvcfs)} \
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
