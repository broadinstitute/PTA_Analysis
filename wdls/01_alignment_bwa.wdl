version 1.0

workflow AlignReadsWorkflow {
    input {
        File read1                              # First-end FASTQ file
        File? read2                            # Second-end FASTQ file (optional)
        File referenceFasta                    # Reference genome FASTA file
        Array[File] bwaIndexFiles              # Array of BWA index files
        String outputPrefix                    # Output file prefix
        String readGroup                       # Read group string

        # Resource settings
        Int threads = 16                       # Number of threads for BWA
        Int memoryGb = 64                      # Memory allocated in GB
        Int diskGb = 200                       # Disk size in GB
        Int preemptible = 3                    # Number of preemptible attempts
        String dockerImage = "us.gcr.io/broad-dsp-lrma/sr-utils:0.2.2"  # Docker image
    }

    call AlignReads {
        input:
            read1 = read1,
            read2 = read2,
            referenceFasta = referenceFasta,
            bwaIndexFiles = bwaIndexFiles,
            outputPrefix = outputPrefix,
            readGroup = readGroup,
            threads = threads,
            memoryGb = memoryGb,
            diskGb = diskGb,
            preemptible = preemptible,
            dockerImage = dockerImage
    }

    # Optional: Only run sorting and fixmate if unsorted BAM exists
    if (defined(AlignReads.outputUnsortedBam)) {
        call PrepareBam {
            input:
                unsortedBam = AlignReads.outputUnsortedBam,
                outputPrefix = outputPrefix,
                threads = threads,
                memoryGb = memoryGb,
                diskGb = diskGb,
                preemptible = preemptible,
                dockerImage = dockerImage
        }
    }

    output {
        File unsortedBam = AlignReads.outputUnsortedBam
        File? sortedBam = PrepareBam.sortedBam
        File? sortedBai = PrepareBam.sortedBai
    }
}

task AlignReads {
    input {
        File read1
        File? read2
        File referenceFasta
        Array[File] bwaIndexFiles
        String outputPrefix
        String readGroup                      # Input the read group string

        # Resource settings
        Int threads
        Int memoryGb
        Int diskGb
        Int preemptible
        String dockerImage
    }

    command <<<
        set -euxo pipefail

        # Align reads and create unsorted BAM file
        bwa mem \
          -K 100000000 \
          -t ~{threads} \
          -Y \
          -R '~{readGroup}' \
          -c 100 \
          -M \
          ~{referenceFasta} \
          ~{read1} \
          ~{if defined(read2) then read2 else ""} \
        | samtools view -b -o ~{outputPrefix}.unsorted.bam -
    >>>

    output {
        File outputUnsortedBam = "~{outputPrefix}.unsorted.bam"
    }

    runtime {
        cpu: threads
        memory: "~{memoryGb} GiB"
        disks: "local-disk ~{diskGb} HDD"
        preemptible: preemptible
        docker: dockerImage
    }
}

task PrepareBam {
    input {
        File unsortedBam
        String outputPrefix

        # Resource settings
        Int threads
        Int memoryGb
        Int diskGb
        Int preemptible
        String dockerImage
    }

    command <<<
        set -euxo pipefail

        # Queryname sort the BAM file
        samtools sort -n -@ ~{threads} \
          -o ~{outputPrefix}.queryname_sorted.bam ~{unsortedBam}

        # Fixmate to correct mate information
        samtools fixmate -m ~{outputPrefix}.queryname_sorted.bam ~{outputPrefix}.fixed.bam

        # Sort BAM file by genomic coordinates
        samtools sort -@ ~{threads} \
          -o ~{outputPrefix}.bam ~{outputPrefix}.fixed.bam

        # Index the sorted BAM file
        samtools index ~{outputPrefix}.bam
    >>>

    output {
        File sortedBam = "~{outputPrefix}.bam"
        File sortedBai = "~{outputPrefix}.bam.bai"
    }

    runtime {
        cpu: threads
        memory: "~{memoryGb} GiB"
        disks: "local-disk ~{diskGb} HDD"
        preemptible: preemptible
        docker: dockerImage
    }
}
