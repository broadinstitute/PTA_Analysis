version 1.0

workflow CreateBwaMem2Index {
    input {
        File reference_fasta  # Input reference genome FASTA file
        String prefix = "reference"  # Prefix for the output index files
        Int? cpu_cores = 16  # Optional: Number of CPU cores for indexing
        Int? memory_gb = 64  # Optional: Memory in GB for the task
        Int? disk_gb = 50  # Optional: Disk space in GB for the task
    }

    call BwaMem2Index {
        input:
            fasta = reference_fasta,
            prefix = prefix,
            cpu_cores = cpu_cores,
            memory_gb = memory_gb,
            disk_gb = disk_gb
    }

    output {
        File amb = BwaMem2Index.amb
        File ann = BwaMem2Index.ann
        File bwt_64 = BwaMem2Index.bwt_64
        File pac = BwaMem2Index.pac
        File sa = BwaMem2Index.sa
        File o123 = BwaMem2Index.o123
    }
}

task BwaMem2Index {
    input {
        File fasta  # Reference FASTA file
        String prefix  # Prefix for output files
        Int? cpu_cores
        Int? memory_gb
        Int? disk_gb
    }

    command <<<

        set -euxo pipefail

        # Ensure the input FASTA file is in the current directory
        ln -s ~{fasta} ~{prefix}.fa

        # Run bwa-mem2 index command
        bwa-mem2 index -p ~{prefix} ~{prefix}.fa
    >>>

    output {
        File amb = "~{prefix}.fa.amb"
        File ann = "~{prefix}.fa.ann"
        File bwt_64 = "~{prefix}.fa.bwt.2bit.64"
        File pac = "~{prefix}.fa.pac"
        File sa = "~{prefix}.fa.sa"
        File o123 = "~{prefix}.fa.0123"
    }

    runtime {
        cpu: select_first([cpu_cores, 8])
        memory: select_first([memory_gb, 64]) + " GiB"
        disks: "local-disk " + select_first([disk_gb, 50]) + " HDD"
        docker: "us.gcr.io/broad-dsp-lrma/sr-utils:0.2.2"
    }
}
