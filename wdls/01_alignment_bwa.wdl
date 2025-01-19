version 1.0

workflow AlignAndIndexWorkflow {

    input {
        File read1                              # First-end FASTQ file
        File? read2                            # Second-end FASTQ file (optional)
        File referenceFasta                    # Reference genome FASTA file
        Array[File] bwaIndexFiles              # Array of BWA index files
        String outputPrefix                    # Output file prefix
        String? readGroup                      # Read group string (optional)
        Int threads = 4                        # Number of threads for BWA
        Int memoryGb = 16                      # Memory allocated in GB
        Int diskGb = 100                       # Disk size in GB
        String dockerImage = "quay.io/biocontainers/bwa:0.7.17--hed695b0_7"  # Docker image
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
            dockerImage = dockerImage
    }

    output {
        File outputBam = AlignReads.outputBam
        File outputBai = AlignReads.outputBai
    }
}

task AlignReads {
    input {
        File read1
        File? read2
        File referenceFasta
        Array[File] bwaIndexFiles
        String outputPrefix
        String? readGroup
        Int threads = 4
        Int memoryGb = 16
        Int diskGb = 100
        String dockerImage = "quay.io/biocontainers/bwa:0.7.17--hed695b0_7"
    }

    command <<<
        set -euxo pipefail
        
        # Run BWA and sort the BAM file
        bwa mem \
          -t ~{threads} \
          ~{if defined(readGroup) then "-R '" + readGroup + "'" else ""} \
          ~{referenceFasta} \
          ~{read1} \
          ~{if defined(read2) then read2 else ""} \
        | samtools sort -@~{threads} -m ~{memoryGb / threads}G -o ~{outputPrefix}.bam -

        # Index the sorted BAM file
        samtools index ~{outputPrefix}.bam
    >>>

    output {
        File outputBam = "~{outputPrefix}.bam"
        File outputBai = "~{outputPrefix}.bam.bai"
    }

    runtime {
        cpu: threads
        memory: "~{memoryGb} GiB"
        disks: "local-disk ~{diskGb} HDD"
        docker: dockerImage
    }

    parameter_meta {
        read1: {description: "First-end FASTQ file"}
        read2: {description: "Second-end FASTQ file (optional)"}
        referenceFasta: {description: "Reference genome FASTA file"}
        bwaIndexFiles: {description: "Array of BWA index files"}
        outputPrefix: {description: "Output file prefix"}
        readGroup: {description: "Read group string (optional)"}
        threads: {description: "Number of threads for BWA"}
        memoryGb: {description: "Memory allocated in GB"}
        diskGb: {description: "Disk size allocated in GB"}
        dockerImage: {description: "Docker image for running BWA"}
    }
}
