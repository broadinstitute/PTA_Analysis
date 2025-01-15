task bwa_mem {
  input {
    File fastq1
    File fastq2
    File reference
    String sample_id
  }

  command {
  }

  output {
    File sam = "~{sample_id}.sam"
  }

  runtime {
    docker: "biocontainers/bwa:v0.7.17"
    cpu: 8
    memory: "32G"
  }
}
task MapReads {
  input {
    File fastq1
    File fastq2
    File reference_fasta
    Int preemptible_tries
    Int threads = 4
    String output_basename
    String docker = "biocontainers/bwa:v0.7.17_cv1"
  }

  command {
    # bwa mem -c 100 -M ~{reference} ~{fastq1} ~{fastq2} > ~{sample_id}.sam
    bwa mem -t ~{threads} ~{reference_fasta} ~{fastq1} ~{fastq2} > ~{output_basename}.sam
    samtools view -bS ~{output_basename}.sam > ~{output_basename}.bam
    samtools sort ~{output_basename}.bam -o ~{output_basename}_sorted.bam
    samtools index ~{output_basename}_sorted.bam
  }

  runtime {
    docker: docker
    memory: "8 GiB"
    cpu: threads
    disks: "local-disk 50 HDD"
  }

  output {
    File sorted_bam = "~{output_basename}_sorted.bam"
    File bam_index = "~{output_basename}_sorted.bam.bai"
  }
}

