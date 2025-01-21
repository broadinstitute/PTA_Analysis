version 1.0

workflow alignment_workflow {
    input {
        File fastq1                              # Input FASTQ file 1
        File fastq2                              # Input FASTQ file 2
        File reference                           # Reference genome file
        String prefix                            # Prefix for output file
        Array[File] bwaIndexes                   # BWA index files
        Int threads = 16                         # Number of threads (default: 16)
        String read_group                        # Read group string

        # Resource settings
        Int memoryGb = 16                        # Memory allocation in GB (default: 16 GB)
        Int diskGb = 200                         # Disk allocation in GB (default: 200 GB)
        Int preemptible = 3                      # Number of preemptible attempts
        String dockerImage = "us.gcr.io/broad-dsp-lrma/sr-utils:0.2.2"  # Docker image
    }

    call bwa_mem_alignment {
        input:
            fastq1 = fastq1,
            fastq2 = fastq2,
            reference = reference,
            prefix = prefix,
            threads = threads,
            read_group = read_group,
            memoryGb = memoryGb,
            diskGb = diskGb,
            preemptible = preemptible,
            dockerImage = dockerImage,
            bwaIndexes = bwaIndexes
    }

    output {
        File bam = bwa_mem_alignment.bam
    }
}

task bwa_mem_alignment {
    input {
        File fastq1                              # Input FASTQ file 1
        File fastq2                              # Input FASTQ file 2
        File reference                           # Reference genome file
        String prefix                            # Prefix for output file
        Int threads                              # Number of threads to use
        String read_group                        # Read group string
        Array[File] bwaIndexes                   # BWA index files

        # Resource settings
        Int memoryGb                             # Memory allocation in GB
        Int diskGb                               # Disk allocation in GB
        Int preemptible                          # Number of preemptible attempts
        String dockerImage                       # Docker image
    }

    command <<<
        set -euxo pipefail

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
    >>>

    output {
        File bam = "~{prefix}_unsorted.bam"
    }

    runtime {
        cpu: threads
        memory: "~{memoryGb} GiB"   # Convert memory to the appropriate format
        disks: "local-disk ~{diskGb} HDD"
        preemptible: preemptible
        docker: dockerImage
    }
}
