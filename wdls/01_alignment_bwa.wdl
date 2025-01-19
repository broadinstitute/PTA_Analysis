version 1.0

workflow SimplifiedAlignWorkflow {
    input {
        File read1                              # First-end FASTQ file
        File? read2                            # Second-end FASTQ file (optional)
        File referenceFasta                    # Reference genome FASTA file
        Array[File] bwaIndexFiles              # Array of BWA index files
        String outputPrefix                    # Output file prefix
        String sampleName                      # Sample name for read group
        String platform = "ILLUMINA"           # Sequencing platform
        String libraryName                     # Library name

        # Resource settings as inputs
        Int threads = 16                       # Number of threads for BWA
        Int memoryGb = 64                      # Memory allocated in GB
        Int diskGb = 200                       # Disk size in GB
        Int preemptible = 0                    # Number of preemptible attempts
        String dockerImage = "us.gcr.io/broad-dsp-lrma/sr-utils:0.2.2"  # Docker image
    }

    call AlignReads {
        input:
            read1 = read1,
            read2 = read2,
            referenceFasta = referenceFasta,
            outputPrefix = outputPrefix,
            sampleName = sampleName,
            platform = platform,
            libraryName = libraryName,
            threads = threads,
            memoryGb = memoryGb,
            diskGb = diskGb,
            preemptible = preemptible,
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
        String outputPrefix
        String sampleName
        String platform
        String libraryName

        # Resource settings
        Int threads
        Int memoryGb
        Int diskGb
        Int preemptible
        String dockerImage
    }

    command <<<
        set -euxo pipefail
        
        # Align reads and create BAM file
        bwa mem \
          -K 100000000 \
          -t ~{threads} \
          -Y \
          -R '@RG\tID:~{sampleName}\tSM:~{sampleName}\tPL:~{platform}\tLB:~{libraryName}' \
          -c 100 \
          -M \
          ~{referenceFasta} \
          ~{read1} \
          ~{if defined(read2) then read2 else ""} \
        | samtools view -b -o ~{outputPrefix}.bam -
        
        # Index the BAM file
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
        preemptible: preemptible
        docker: dockerImage
    }

    parameter_meta {
        threads: {description: "Number of threads for BWA"}
        memoryGb: {description: "Memory allocated in GB"}
        diskGb: {description: "Disk size allocated in GB"}
        preemptible: {description: "Number of preemptible attempts"}
        dockerImage: {description: "Docker image for running BWA"}
    }
}
