version 1.0

task Index {
    input {
        String donor_id
        String sample_id
        File bam_file
        Int threads = 1
        
        # Runtime parameters
        String docker = "quay.io/biocontainers/sambamba:0.8.2--h98b6b92_2"
        Int memory_gb = 4
        Int disk_size_gb = 100
    }

    command <<<
        set -euo pipefail
        sambamba index -t ~{threads} ~{bam_file} ~{bam_file}.bai
    >>>

    output {
        File bai_file = "~{bam_file}.bai"
    }

    runtime {
        docker: docker
        memory: "~{memory_gb}GB"
        cpu: threads
        disks: "local-disk ~{disk_size_gb} SSD"
    }

    meta {
        author: "Converted from Nextflow"
        description: "Creates an index file for a BAM file using Sambamba"
    }
}

workflow SamtoolsIndex {
    input {
        String donor_id
        String sample_id
        File bam_file
    }

    call Index {
        input:
            donor_id = donor_id,
            sample_id = sample_id,
            bam_file = bam_file
    }

    output {
        File bai_file = Index.bai_file
    }
}