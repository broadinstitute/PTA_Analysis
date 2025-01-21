version 1.0

workflow SortBamWorkflow {
    input {
        File bamFile                     # Input BAM file to be sorted
        String sortType = "coordinate"  # Type of sorting: "coordinate" or "queryname"
        Int threads = 4                 # Number of threads to use
        Int memoryGb = 16               # Memory allocation in GB
        Int diskGb = 100                # Disk size in GB
        String dockerImage = "us.gcr.io/broad-dsp-lrma/sr-utils:0.2.2" # Docker image
        String sampleName               # Name of the sample (used in output file)
    }

    call SortBam {
        input:
            bamFile = bamFile,
            sortType = sortType,
            threads = threads,
            memoryGb = memoryGb,
            diskGb = diskGb,
            dockerImage = dockerImage,
            sampleName = sampleName
    }

    output {
        File sortedBam = SortBam.sortedBam
        File? sortedBai = SortBam.sortedBai
    }
}

task SortBam {
    input {
        File bamFile                     # Input BAM file
        String sortType                  # Type of sorting: "coordinate" or "queryname"
        Int threads                      # Number of threads
        Int memoryGb                     # Memory allocation in GB
        Int diskGb                       # Disk size in GB
        String dockerImage               # Docker image
        String sampleName                # Name of the sample (used in output file)
    }

    command <<<
        set -euxo pipefail

        echo "Using Docker image: ~{dockerImage}"
        echo "Sorting BAM file: ~{bamFile} with sort type: ~{sortType}"
        
        # Print samtools version for debugging
        samtools --version

        # Determine sorting command based on user input
        if [[ "~{sortType}" == "coordinate" ]]; then
            echo "Performing coordinate sorting..."
            samtools sort -@ ~{threads} \
              -o ~{sampleName}.sorted.bam ~{bamFile}
        elif [[ "~{sortType}" == "queryname" ]]; then
            echo "Performing queryname sorting..."
            samtools sort -n -@ ~{threads} \
              -o ~{sampleName}.sorted.bam ~{bamFile}
        else
            echo "ERROR: Invalid sort type: ~{sortType}. Use 'coordinate' or 'queryname'."
            exit 1
        fi

        # Index the BAM file if coordinate sorting is performed
        if [[ "~{sortType}" == "coordinate" ]]; then
            echo "Indexing coordinate-sorted BAM..."
            samtools index ~{sampleName}.sorted.bam
        fi
    >>>

    output {
        File sortedBam = "~{sampleName}.sorted.bam"
        File? sortedBai = "~{sampleName}.sorted.bam.bai"  # Only produced if coordinate sorting
    }

    runtime {
        cpu: threads
        memory: "~{memoryGb} GiB"
        disks: "local-disk ~{diskGb} HDD"
        docker: dockerImage
    }
}
