version 1.0

workflow VariantCalling {
  input {
    File reference_fasta
    File reference_fasta_index
    File input_bam
    File input_bam_index
    String output_vcf_basename
    Int preemptible_tries = 1
    Int memory_gb = 8
    Int disk_gb = 50
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.3.0.0"
  }

  call HaplotypeCaller {
    input:
      reference_fasta = reference_fasta,
      reference_fasta_index = reference_fasta_index,
      input_bam = input_bam,
      input_bam_index = input_bam_index,
      output_vcf_basename = output_vcf_basename,
      preemptible_tries = preemptible_tries,
      memory_gb = memory_gb,
      disk_gb = disk_gb,
      gatk_docker = gatk_docker
  }

  output {
    File raw_variants_vcf = HaplotypeCaller.raw_variants_vcf
    File raw_variants_vcf_index = HaplotypeCaller.raw_variants_vcf_index
  }
}

### Task for HaplotypeCaller
task HaplotypeCaller {
  input {
    File reference_fasta
    File reference_fasta_index
    File input_bam
    File input_bam_index
    String output_vcf_basename
    Int preemptible_tries
    Int memory_gb
    Int disk_gb
    String gatk_docker
  }

  command {
    gatk --java-options "-Xmx~{memory_gb}G" HaplotypeCaller \
      -R ~{reference_fasta} \
      -I ~{input_bam} \
      -O ~{output_vcf_basename}.vcf \
      --emit-ref-confidence GVCF
  }

  runtime {
    docker: gatk_docker
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: preemptible_tries
  }

  output {
    File raw_variants_vcf = "~{output_vcf_basename}.vcf"
    File raw_variants_vcf_index = "~{output_vcf_basename}.vcf.idx"
  }
}