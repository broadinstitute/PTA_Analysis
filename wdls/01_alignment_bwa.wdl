version 1.0

workflow AlignAndIndexWorkflow {
    input {
        File read1                              # First-end FASTQ file
        File? read2                            # Second-end FASTQ file (optional)
        File referenceFasta                    # Reference genome FASTA file
        Array[File] bwaIndexFiles              # Array of BWA index files
        String outputPrefix                    # Output file prefix
        String sample_name                     # Sample name for read group
        String? readGroup                      # Read group string (optional)

        # Resource settings as inputs
        Int threads                             # Number of threads for BWA
        Int memoryGb                            # Memory allocated in GB
        Int diskGb                              # Disk size in GB
        String dockerImage                      # Docker image to use
    }

    call ExtractReadGroup {
        input:
            fastq = read1,
            sample_name = sample_name
    }

    call AlignReads {
        input:
            read1 = read1,
            read2 = read2,
            referenceFasta = referenceFasta,
            bwaIndexFiles = bwaIndexFiles,
            outputPrefix = outputPrefix,
            readGroup = ExtractReadGroup.read_group,
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

task ExtractReadGroup {
    input {
        File fastq
        String sample_name
    }

    command <<<
        zcat ~{fastq} | head -n 1 | sed 's/^@//' | awk -F':' '{printf "@RG\\tID:%s_%s\\tPL:illumina\\tLB:%s\\tSM:%s\\n", $1, $3, $NF, "~{sample_name}"}'
    >>>

    output {
        String read_group = read_string(stdout())
    }

    runtime {
        cpu: 1
        memory: "1G"
        disk: "10G"
        docker: "us.gcr.io/broad-dsp-lrma/sr-utils:0.2.2"
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

        # Resource settings
        Int threads
        Int memoryGb
        Int diskGb
        String dockerImage
    }

    command <<<
        set -euxo pipefail
        
        bwa mem \
          -t ~{threads} \
          ~{if defined(readGroup) then "-R '" + readGroup + "'" else ""} \
          ~{referenceFasta} \
          ~{read1} \
          ~{if defined(read2) then read2 else ""} \
        | samtools sort -@~{threads} -m ~{memoryGb / threads}G -o ~{outputPrefix}.bam -

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
        threads: {description: "Number of threads for BWA"}
        memoryGb: {description: "Memory allocated in GB"}
        diskGb: {description: "Disk size allocated in GB"}
        dockerImage: {description: "Docker image for running BWA"}
    }
}
