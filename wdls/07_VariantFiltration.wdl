version 1.0

# Workflow: VariantFiltration
# Description: This workflow filters raw variant calls from GATK HaplotypeCaller (with EMIT_ALL_CONFIDENT_SITES)
# Author: Constantijn Scharlee
# Date: 2025-02-23
# using GATK VariantFiltration, following the filtering criteria outlined in the reference paper.

workflow VariantFiltration {
  parameter_meta {
    reference_fasta: "Path to the reference genome in FASTA format."
    reference_fasta_index: "Path to the index file for the reference FASTA."
    reference_dict: "Path to the sequence dictionary for the reference FASTA."
    input_vcf: "Path to the input VCF file (must be from HaplotypeCaller with EMIT_ALL_CONFIDENT_SITES)."
    input_vcf_index: "Path to the index file for the input VCF."
    output_vcf_basename: "Base name for the output VCF file (compressed)."
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
    File input_vcf
    File input_vcf_index
    String output_vcf_basename
    Int preemptible_tries = 1
    Int memory_gb = 32
    Int disk_gb = 200
    Int cpu = 4
    String gatk_docker = "broadinstitute/gatk:4.6.1.0"
  }

  call VariantFiltrationTask {
    input:
      reference_fasta = reference_fasta,
      reference_fasta_index = reference_fasta_index,
      reference_dict = reference_dict,
      input_vcf = input_vcf,
      input_vcf_index = input_vcf_index,
      output_vcf_basename = output_vcf_basename,
      preemptible_tries = preemptible_tries,
      memory_gb = memory_gb,
      disk_gb = disk_gb,
      cpu = cpu,
      gatk_docker = gatk_docker
  }

  output {
    File filtered_vcf = VariantFiltrationTask.filtered_vcf
    File filtered_vcf_index = VariantFiltrationTask.filtered_vcf_index
  }
}

### **Task for VariantFiltration**
task VariantFiltrationTask {
  parameter_meta {
    reference_fasta: "Path to the reference genome in FASTA format."
    reference_fasta_index: "Path to the index file for the reference FASTA."
    reference_dict: "Path to the sequence dictionary for the reference FASTA."
    input_vcf: "Path to the input VCF file (must be from HaplotypeCaller with EMIT_ALL_CONFIDENT_SITES)."
    input_vcf_index: "Path to the index file for the input VCF."
    output_vcf_basename: "Base name for the output VCF file (compressed)."
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
    File input_vcf
    File input_vcf_index
    String output_vcf_basename
    Int preemptible_tries
    Int memory_gb
    Int disk_gb
    Int cpu
    String gatk_docker
  }

  command {
    gatk --java-options "-Xmx~{memory_gb}G" VariantFiltration \
      -R ~{reference_fasta} \
      -V ~{input_vcf} \
      -O ~{output_vcf_basename}.vcf.gz \
      --filter-expression "QD < 2.0" --filter-name "SNP_LowQualityDepth" \
      --filter-expression "MQ < 40.0" --filter-name "SNP_MappingQuality" \
      --filter-expression "FS > 60.0" --filter-name "SNP_StrandBias" \
      --filter-expression "HaplotypeScore > 13.0" --filter-name "SNP_HaplotypeScoreHigh" \
      --filter-expression "MQRankSum < -12.5" --filter-name "SNP_MQRankSumLow" \
      --filter-expression "ReadPosRankSum < -8.0" --filter-name "SNP_ReadPosRankSumLow" \
      --filter-expression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" --filter-name "SNP_HardToValidate" \
      --filter-expression "DP < 5" --filter-name "SNP_LowCoverage" \
      --filter-expression "QUAL < 30" --filter-name "SNP_VeryLowQual" \
      --filter-expression "QUAL >= 30.0 && QUAL < 50.0" --filter-name "SNP_LowQual" \
      --filter-expression "SOR > 4.0" --filter-name "SNP_SOR" \
      -cluster 3 -window 10
  }

  runtime {
    docker: gatk_docker
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: preemptible_tries
    cpu: cpu
  }

  output {
    File filtered_vcf = "~{output_vcf_basename}.vcf.gz"
    File filtered_vcf_index = "~{output_vcf_basename}.vcf.gz.tbi"
  }
}

### **Broad Institute License**
# Copyright (c) 2025 Broad Institute
# This software is distributed under the BSD-3-Clause License.
