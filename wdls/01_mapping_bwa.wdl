version 1.0

task bwa_mem_alignment {
    input {
        File fastq1               # Input FASTQ file 1
        File fastq2               # Input FASTQ file 2
        File reference            # Reference genome file
        String prefix             # Prefix for output file
        Int threads               # Number of threads to use
        String read_group         # Read group string
        String memory             # Memory allocation (e.g., "16 GB")
        String disk               # Disk allocation (e.g., "200 GB")
    }

    command {
        echo "Running bwa mem..."
        bwa mem \
            -K 100000000 \
            -t ~{threads} \
            -Y \
            -R "~{read_group}" \
            -c 100 \
            -M \
            ~{reference} \
            ~{fastq1} ~{fastq2} | \
            samtools view -b -o ~{prefix}_unsorted.bam -
    }

    output {
        File bam = "~{prefix}_unsorted.bam"
    }

    runtime {
        cpu: threads
        memory: memory
        disks: "local-disk ~{disk}"#
        docker: "us.gcr.io/broad-dsp-lrma/sr-utils:0.2.2"
    }
}

workflow alignment_workflow {
    input {
        File fastq1               # Input FASTQ file 1
        File fastq2               # Input FASTQ file 2
        File reference            # Reference genome file
        String prefix             # Prefix for output file
        Int threads = 16          # Number of threads (default: 16)
        String read_group         # Read group string
        String memory = "16 GB"   # Memory allocation (default: 16 GB)
        String disk = "200 GB"    # Disk allocation (default: 200 GB)
    }

    call bwa_mem_alignment {
        input:
            fastq1 = fastq1,
            fastq2 = fastq2,
            reference = reference,
            prefix = prefix,
            threads = threads,
            read_group = read_group,
            memory = memory,
            disk = disk
    }

    output {
        File bam = bwa_mem_alignment.bam
    }
}
