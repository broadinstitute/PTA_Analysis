version 1.0

# Workflow: VariantCalling
# Author: Shadi Zaheri
# Date: 2025-01-26
# Description: This workflow performs single-sample variant calling on whole genome sequencing (WGS) data 
# using GATK HaplotypeCaller, following the settings from the reference paper.
# The input BAM should be the output from Realignment, and the workflow assumes known variant databases are available.

workflow VariantCalling {
  parameter_meta {
    reference_fasta: "Path to the reference genome in FASTA format."
    reference_fasta_index: "Path to the index file for the reference FASTA."
    reference_dict: "Path to the sequence dictionary for the reference FASTA."
    input_bam: "Path to the input BAM file."
    input_bam_index: "Path to the index file for the input BAM."
    output_vcf_basename: "Base name for the output VCF file (compressed)."
    interval_list: "Optional: Path to an interval list file (e.g., BED or Picard format) for targeted variant calling."
    preemptible_tries: "Number of preemptible retries allowed for each task."
    memory_gb: "Memory allocated for each task in gigabytes."
    disk_gb: "Disk space allocated for each task in gigabytes."
    cpu: "Number of CPU cores allocated for each task."
    gatk_docker: "Docker image for GATK tools (default: broadinstitute/gatk:4.6.1.0)."
  }

  input {
    File reference_fasta
    File reference_fasta_index
    File reference_dict  # Dictionary file (.dict)
    File input_bam
    File input_bam_index
    String output_vcf_basename
    File? interval_list  # Optional input
    Int preemptible_tries = 1
    Int memory_gb = 32  # Increased for large WGS datasets
    Int disk_gb = 200  # Increased to handle large files
    Int cpu = 4
    String gatk_docker = "broadinstitute/gatk:4.6.1.0"  # Matching the paper
  }

  call HaplotypeCaller {
    input:
      reference_fasta = reference_fasta,
      reference_fasta_index = reference_fasta_index,
      reference_dict = reference_dict,
      input_bam = input_bam,
      input_bam_index = input_bam_index,
      output_vcf_basename = output_vcf_basename,
      interval_list = interval_list,  # Optional interval list
      preemptible_tries = preemptible_tries,
      memory_gb = memory_gb,
      disk_gb = disk_gb,
      cpu = cpu,
      gatk_docker = gatk_docker
  }

  output {
    File raw_variants_vcf = HaplotypeCaller.raw_variants_vcf
    File raw_variants_vcf_index = HaplotypeCaller.raw_variants_vcf_index
  }
}

### **Task for HaplotypeCaller**
task HaplotypeCaller {
  parameter_meta {
    reference_fasta: "Path to the reference genome in FASTA format."
    reference_fasta_index: "Path to the index file for the reference FASTA."
    reference_dict: "Path to the sequence dictionary for the reference FASTA."
    input_bam: "Path to the input BAM file."
    input_bam_index: "Path to the index file for the input BAM."
    output_vcf_basename: "Base name for the output VCF file (compressed)."
    interval_list: "Optional: Path to an interval list file (e.g., BED or Picard format) for targeted variant calling."
    preemptible_tries: "Number of preemptible retries allowed for each task."
    memory_gb: "Memory allocated for each task in gigabytes."
    disk_gb: "Disk space allocated for each task in gigabytes."
    cpu: "Number of CPU cores allocated for each task."
    gatk_docker: "Docker image for GATK tools (default: broadinstitute/gatk:4.6.1.0)."
  }

  input {
    File reference_fasta
    File reference_fasta_index
    File input_bam
    File input_bam_index
    File reference_dict  
    String output_vcf_basename
    File? interval_list  # Optional input
    Int preemptible_tries
    Int memory_gb
    Int disk_gb
    Int cpu
    String gatk_docker
  }

  command {
    gatk --java-options "-Xmx~{memory_gb}G" HaplotypeCaller \
      -R ~{reference_fasta} \
      -I ~{input_bam} \
      -O ~{output_vcf_basename}.vcf.gz \
      --emit-ref-confidence GVCF \
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
    File raw_variants_vcf = "~{output_vcf_basename}.vcf.gz"
    File raw_variants_vcf_index = "~{output_vcf_basename}.vcf.gz.tbi"
  }
}
