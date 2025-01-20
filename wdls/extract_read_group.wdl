version 1.0

workflow ExtractReadGroupWorkflow {
    input {
        File read1                              # First-end FASTQ file
        String sampleName                      # Sample name for the read group
        String platform = "ILLUMINA"           # Sequencing platform
        String? libraryName                    # Optional library name

        # Runtime attributes
        Int cpu = 1
        Int memoryGb = 4
        Int diskGb = 20
        Int preemptible = 0
        String dockerImage = "us.gcr.io/broad-dsp-lrma/sr-utils:0.2.2"
    }

    call ExtractReadGroup {
        input:
            fastq = read1,
            sampleName = sampleName,
            platform = platform,
            libraryName = libraryName,
            cpu = cpu,
            memoryGb = memoryGb,
            diskGb = diskGb,
            preemptible = preemptible,
            dockerImage = dockerImage
    }

    output {
        String readGroup = ExtractReadGroup.read_group
    }
}

task ExtractReadGroup {
    input {
        File fastq
        String sample_name
    }

    command <<<
        zcat ~{fastq} | head -n 1 | awk -F':' '{printf "@RG\\tID:%s_%s\\tPL:illumina\\tLB:%s\\tSM:%s", $1, $3, $NF, "~{sample_name}"}'
    >>>

    output {
        String read_group = read_string(stdout())
    }

    runtime {
        cpu: 1
        memory: "4G"
        disk: "10G"
        docker: "us.gcr.io/broad-dsp-lrma/sr-utils:0.2.2"
    }
}
