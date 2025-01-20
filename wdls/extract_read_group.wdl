version 1.0

workflow ExtractReadGroupWorkflow {
    input {
        File read1                              # First-end FASTQ file
        String sample_name                      # Sample name for the read group
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
            sample_name = sample_name,
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
        String platform
        String? libraryName

        # Runtime attributes
        Int cpu
        Int memoryGb
        Int diskGb
        Int preemptible
        String dockerImage
    }

    command <<<
        zcat ~{fastq} | head -n 1 | awk -F':' '{printf "@RG\\tID:%s_%s\\tPL:illumina\\tLB:%s\\tSM:%s", $1, $3, $NF, "~{sample_name}"}'
    >>>

    output {
        String read_group = read_string(stdout())
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGb} GiB"
        disks: "local-disk ~{diskGb} HDD"
        preemptible: preemptible
        docker: dockerImage
    }
}
