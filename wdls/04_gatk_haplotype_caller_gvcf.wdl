version 1.0

# Workflow: GVCFCalling
# Author: Shadi Zaheri
# Date: 2025-05-01
# Description: This workflow performs per-sample variant calling using GATK HaplotypeCaller in GVCF mode.
# It generates a .g.vcf.gz file per sample to enable scalable joint genotyping across multiple samples.

workflow GVCFCalling {
  parameter_meta {
    reference_fasta: "Path to the reference genome in FASTA format."
    reference_fasta_index: "Path to the index file for the reference FASTA."
    reference_dict: "Path to the sequence dictionary for the reference FASTA."
    input_bam: "Path to the input BAM file."
    input_bam_index: "Path to the index file for the input BAM."
    output_gvcf_basename: "Base name for the output GVCF file (compressed)."
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
    File reference_dict  
    File input_bam
    File input_bam_index
    String output_gvcf_basename
    File? interval_list  
    Int preemptible_tries = 1
    Int memory_gb = 32 
    Int disk_gb = 200  
    Int cpu = 4
    String gatk_docker = "broadinstitute/gatk:4.6.1.0"  
  }

  call HaplotypeCallerGVCF {
    input:
      reference_fasta = reference_fasta,
      reference_fasta_index = reference_fasta_index,
      reference_dict = reference_dict,
      input_bam = input_bam,
      input_bam_index = input_bam_index,
      output_gvcf_basename = output_gvcf_basename,
      interval_list = interval_list,
      preemptible_tries = preemptible_tries,
      memory_gb = memory_gb,
      disk_gb = disk_gb,
      cpu = cpu,
      gatk_docker = gatk_docker
  }

  output {
    File gvcf = HaplotypeCallerGVCF.gvcf
    File gvcf_index = HaplotypeCallerGVCF.gvcf_index
  }
}

### Task: HaplotypeCallerGVCF
task HaplotypeCallerGVCF {
  parameter_meta {
    reference_fasta: "Path to the reference genome in FASTA format."
    reference_fasta_index: "Path to the index file for the reference FASTA."
    reference_dict: "Path to the sequence dictionary for the reference FASTA."
    input_bam: "Path to the input BAM file."
    input_bam_index: "Path to the index file for the input BAM."
    output_gvcf_basename: "Base name for the output GVCF file (compressed)."
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
    File reference_dict
    File input_bam
    File input_bam_index
    String output_gvcf_basename
    File? interval_list
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
      -O ~{output_gvcf_basename}.g.vcf.gz \
      -ERC GVCF \
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
    File gvcf = "~{output_gvcf_basename}.g.vcf.gz"
    File gvcf_index = "~{output_gvcf_basename}.g.vcf.gz.tbi"
  }
}

