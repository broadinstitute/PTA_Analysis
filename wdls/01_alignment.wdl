version 1.0

task cram_to_fastq {
    input {
        File input_cram
        File reference_fasta
        File reference_fasta_index
        String output_prefix
        Int threads = 4
        Int disk_size = ceil(size(input_cram, "GB") * 3 + 20)
    }

    command <<<
        set -euo pipefail
        
        # Convert CRAM to FASTQ (paired-end)
        samtools fastq \
            -@ ~{threads} \
            -1 ~{output_prefix}_R1.fastq.gz \
            -2 ~{output_prefix}_R2.fastq.gz \
            -0 /dev/null \
            -s /dev/null \
            -T RX \
            -n \
            --reference ~{reference_fasta} \
            ~{input_cram}
    >>>

    output {
        File fastq1 = "~{output_prefix}_R1.fastq.gz"
        File fastq2 = "~{output_prefix}_R2.fastq.gz"
    }

    runtime {
        docker: "staphb/samtools:latest"
        memory: "8G"
        cpu: threads
        disks: "local-disk " + disk_size + " SSD"
        preemptible: 3
    }
}

task bwa_mem {
  input {
    File fastq1
    File fastq2
    File reference_fasta
    File reference_fasta_index
    File reference_fasta_amb
    File reference_fasta_ann
    File reference_fasta_bwt
    File reference_fasta_pac
    File reference_fasta_sa
    String output_prefix
    # Runtime parameters
    Int threads = 8
    String memory = "16G"
    Int disk_size = ceil(size(fastq1, "GB") + size(fastq2, "GB") + size(reference_fasta, "GB") + 20)
  }

  command <<<
    set -euo pipefail
    
    # Run BWA MEM and pipe directly to SAMtools to create sorted BAM
    bwa mem \
      -t ~{threads} \
      -M \
      -R "@RG\tID:~{output_prefix}\tSM:~{output_prefix}\tPL:illumina" \
      ~{reference_fasta} \
      ~{fastq1} ~{fastq2} \
    | samtools sort -@ ~{threads} -O bam -o ~{output_prefix}.bam
    
    # Index the BAM file
    samtools index ~{output_prefix}.bam
  >>>

  output {
    File output_bam = "~{output_prefix}.bam"
    File output_bam_index = "~{output_prefix}.bam.bai"
  }

  runtime {
    docker: "quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:219b6c272b25e7e642ae3ff0bf0c5c81a5135ab4-0" # Contains both BWA and samtools
    memory: memory
    cpu: threads
    disks: "local-disk " + disk_size + " SSD"
    preemptible: 3
    maxRetries: 1
  }
}

workflow bwa_mem_workflow {
    input {
        File input_cram
        File reference_fasta
        File reference_fasta_index
        File reference_fasta_amb
        File reference_fasta_ann
        File reference_fasta_bwt
        File reference_fasta_pac
        File reference_fasta_sa
        String output_prefix
    }

    call cram_to_fastq {
        input:
            input_cram = input_cram,
            reference_fasta = reference_fasta,
            reference_fasta_index = reference_fasta_index,
            output_prefix = output_prefix
    }

    call bwa_mem {
        input:
            fastq1 = cram_to_fastq.fastq1,
            fastq2 = cram_to_fastq.fastq2,
            reference_fasta = reference_fasta,
            reference_fasta_index = reference_fasta_index,
            reference_fasta_amb = reference_fasta_amb,
            reference_fasta_ann = reference_fasta_ann,
            reference_fasta_bwt = reference_fasta_bwt,
            reference_fasta_pac = reference_fasta_pac,
            reference_fasta_sa = reference_fasta_sa,
            output_prefix = output_prefix
    }

    output {
        File aligned_bam = bwa_mem.output_bam
        File aligned_bam_index = bwa_mem.output_bam_index
    }
}
