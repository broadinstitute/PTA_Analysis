version 1.0

workflow CreateBwaMem2IndexWorkflow {
    input {
        File reference_fasta  # Input reference genome FASTA file
        Int? cpu_cores = 16  # Optional: Number of CPU cores for indexing
        Int? memory_gb = 128  # Optional: Memory in GB for the task
        Int? disk_gb = 500  # Optional: Disk space in GB for the task
    }

    call BwaMem2Index {
        input:
            fasta = reference_fasta,
            cpu_cores = cpu_cores,
            memory_gb = memory_gb,
            disk_gb = disk_gb
    }

    output {
        File? amb = BwaMem2Index.amb
        File? ann = BwaMem2Index.ann
        File? bwt_64 = BwaMem2Index.bwt_64
        File? pac = BwaMem2Index.pac
        File? sa = BwaMem2Index.sa
        File? o123 = BwaMem2Index.o123
    }
}

task BwaMem2Index {
    input {
        File fasta  # Reference FASTA file
        Int? cpu_cores
        Int? memory_gb
        Int? disk_gb
    }

    command <<<
        set -euxo pipefail

        # Ensure the input FASTA file is in the current directory
        fasta_basename=$(basename ~{fasta})
        ln -s ~{fasta} ${fasta_basename}

        # Run bwa-mem2 index command
        bwa-mem2 index ${fasta_basename}
    >>>

    output {
        File? amb = "~{basename(fasta)}.amb"
        File? ann = "~{basename(fasta)}.ann"
        File? bwt_64 = "~{basename(fasta)}.bwt.2bit.64"
        File? pac = "~{basename(fasta)}.pac"
        File? sa = "~{basename(fasta)}.sa"
        File? o123 = "~{basename(fasta)}.0123"
    }

    runtime {
        cpu: select_first([cpu_cores, 16])
        memory: select_first([memory_gb, 128]) + " GiB"
        disks: "local-disk " + select_first([disk_gb, 500]) + " HDD"
        docker: "us.gcr.io/broad-dsp-lrma/sr-utils:0.2.2"
    }
}
