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
        String sampleName
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
        zcat ~{fastq} | head -n 1 | sed 's/^@//' | awk -F':' '{printf "@RG\\tID:%s_%s\\tPL:~{platform}\\tSM:%s~{if defined(libraryName) then "\\tLB:" + libraryName else ""}\\n", $1, $3, "~{sampleName}"}'
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
