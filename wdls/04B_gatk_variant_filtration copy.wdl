version 1.0

workflow VariantFiltration {
  input {
    File raw_variants_vcf
    File raw_variants_vcf_index
    File reference_fasta
    File reference_fasta_index
    String output_vcf_basename
    Int preemptible_tries = 1
    Int memory_gb = 4
    Int disk_gb = 20
    Int cpu = 4
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.3.0.0"
  }

  call ApplyVariantFiltration {
    input:
      raw_variants_vcf = raw_variants_vcf,
      raw_variants_vcf_index = raw_variants_vcf_index,
      reference_fasta = reference_fasta,
      reference_fasta_index = reference_fasta_index,
      output_vcf_basename = output_vcf_basename,
      preemptible_tries = preemptible_tries,
      memory_gb = memory_gb,
      disk_gb = disk_gb,
      cpu = cpu,
      gatk_docker = gatk_docker
  }

  output {
    File filtered_variants_vcf = ApplyVariantFiltration.filtered_variants_vcf
    File filtered_variants_vcf_index = ApplyVariantFiltration.filtered_variants_vcf_index
  }
}

task ApplyVariantFiltration {
  input {
    File raw_variants_vcf
    File raw_variants_vcf_index
    File reference_fasta
    File reference_fasta_index
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
      -V ~{raw_variants_vcf} \
      -O ~{output_vcf_basename}.vcf \
      --filter-expression "QD < 2.0" --filter-name "LowQualityDepth" \
      --filter-expression "MQ < 40.0" --filter-name "MappingQuality" \
      --filter-expression "FS > 60.0" --filter-name "StrandBias" \
      --filter-expression "HaplotypeScore > 13.0" --filter-name "HaplotypeScoreHigh" \
      --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSumLow" \
      --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSumLow" \
      --filter-expression "MQ0 >= 4 && ((MQ0/(1.0 * DP)) > 0.1)" --filter-name "HardToValidate" \
      --filter-expression "DP < 5" --filter-name "LowCoverage" \
      --filter-expression "QUAL < 30.0" --filter-name "VeryLowQual" \
      --filter-expression "QUAL >= 30.0 && QUAL < 50.0" --filter-name "LowQual" \
      --filter-expression "SOR > 4.0" --filter-name "StrandOddsRatio"
  }

  runtime {
    docker: gatk_docker
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: preemptible_tries
    cpu: cpu
  }

  output {
    File filtered_variants_vcf = "~{output_vcf_basename}.vcf"
    File filtered_variants_vcf_index = "~{output_vcf_basename}.vcf.idx"
  }
}