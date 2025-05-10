version 1.0

# Copyright (c) 2025 The Broad Institute
# Distributed under the BSD-3-Clause License

workflow VariantFiltration {
  meta {
    author: "Shadi Zaheri"
    description: "This workflow filters raw variant calls produced by GATK GenotypeGVCFs, applying recommended quality thresholds."
  }

  parameter_meta {
    reference_fasta: "Reference genome FASTA file used for alignment and variant calling."
    reference_fasta_index: "Index file for the reference FASTA (.fai)."
    reference_dict: "Sequence dictionary file for the reference genome (.dict)."
    input_vcf: "Input VCF file from joint genotyping step (GenotypeGVCFs)."
    input_vcf_index: "Index file (.tbi) for the input VCF."
    output_vcf_basename: "Base name for the filtered output VCF file (will be compressed and indexed)."
    preemptible_tries: "Number of preemptible retries allowed."
    memory_gb: "Memory allocated to the task in GB."
    disk_gb: "Disk size in GB for local scratch space."
    cpu: "Number of CPUs to allocate to the task."
    gatk_docker: "Docker image used for running GATK (default: broadinstitute/gatk:4.6.1.0)."
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

task VariantFiltrationTask {
  meta {
    author: "Constantijn Scharlee"
    description: "Filters raw variants using GATK VariantFiltration based on quality metrics (QD, MQ, FS, etc.)."
  }

  parameter_meta {
    reference_fasta: "Reference genome FASTA file."
    reference_fasta_index: "Index (.fai) for the reference FASTA."
    reference_dict: "Sequence dictionary for the reference genome."
    input_vcf: "Input joint-genotyped VCF file."
    input_vcf_index: "Index (.tbi) for the input VCF."
    output_vcf_basename: "Base name for filtered output VCF."
    preemptible_tries: "Retry attempts if using preemptible VMs."
    memory_gb: "Memory in GB."
    disk_gb: "Disk size in GB."
    cpu: "CPU cores to use."
    gatk_docker: "Docker image for GATK tool."
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
      --filter-expression "QD < 2.0" --filter-name "LowQD" \
      --filter-expression "MQ < 40.0" --filter-name "LowMQ" \
      --filter-expression "FS > 60.0" --filter-name "HighStrandBias" \
      --filter-expression "MQRankSum < -12.5" --filter-name "LowMQRankSum" \
      --filter-expression "ReadPosRankSum < -8.0" --filter-name "LowReadPosRankSum" \
      --filter-expression "DP < 5" --filter-name "LowCoverage" \
      --filter-expression "QUAL < 30.0" --filter-name "VeryLowQual" \
      --filter-expression "QUAL >= 30.0 && QUAL < 50.0" --filter-name "LowQual" \
      --filter-expression "SOR > 4.0" --filter-name "HighSOR" \
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
