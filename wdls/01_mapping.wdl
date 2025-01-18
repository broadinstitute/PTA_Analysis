version 1.0

workflow BwaMem2Alignment {

    input {
        File fq1           # FASTQ file for read 1
        File fq2           # FASTQ file for read 2
        File ref_fasta     # Reference genome FASTA file
        String read_group  # Read group information (e.g., '@RG\tID:foo\tSM:bar')
        String prefix = "aligned"  # Output file prefix

        Int? cpu_cores     # Number of CPU cores
        Int? memory_gb     # Amount of memory in GB
        Int? disk_gb       # Amount of disk space in GB
    }

    call BwaMem2Task {
        input:
            fq1 = fq1,
            fq2 = fq2,
            ref_fasta = ref_fasta,
            read_group = read_group,
            prefix = prefix,
            cpu_cores = cpu_cores,
            memory_gb = memory_gb,
            disk_gb = disk_gb
    }

    output {
        File aligned_bam = BwaMem2Task.bam
    }
}


task BwaMem2Task {
    input {
        File fq1           # FASTQ file for read 1
        File fq2           # FASTQ file for read 2
        File ref_fasta     # Reference genome FASTA file
        String read_group  # Read group information
        String prefix = "aligned"  # Output file prefix

        Int? cpu_cores     # Number of CPU cores
        Int? memory_gb     # Amount of memory in GB
        Int? disk_gb       # Amount of disk space in GB
    }

    Int default_cpu = 4
    Int default_memory = 16
    Int default_disk = 100

    Int effective_cpu = select_first([cpu_cores, default_cpu])
    Int effective_memory = select_first([memory_gb, default_memory])
    Int effective_disk = select_first([disk_gb, default_disk])

    command <<< 
        set -euxo pipefail

        bwa-mem2 mem \
            -K 100000000 \
            -t ~{effective_cpu} \
            -R '~{read_group}' \
            -c 100 -M \
            ~{ref_fasta} \
            ~{fq1} ~{fq2} | \
        samtools view -bS -o ~{prefix}.bam
    >>>

    output {
        File bam = "~{prefix}.bam"
    }

    runtime {
        cpu: effective_cpu
        memory: effective_memory + " GiB"
        disks: "local-disk " + effective_disk + " HDD"
        docker: "us.gcr.io/broad-dsp-lrma/sr-utils:0.2.2"
    }
}
