version 1.0

task SambambaIndex {
    input {
        String donor_id
        String sample_id
        File bam_file
        Int threads = 4
        String docker = "quay.io/biocontainers/sambamba:0.8.2--h98b6b92_2"
    }

    command <<<
        set -euo pipefail
        sambamba index -t ~{threads} ~{bam_file} ~{bam_file}.bai
    >>>

    runtime {
        docker: docker
        cpu: threads
        memory: "4 GB"
    }

    output {
        File bai_file = "~{bam_file}.bai"
    }
}

workflow Index {
    input {
        String donor_id
        String sample_id
        File bam_file
    }

    call SambambaIndex {
        input:
            donor_id = donor_id,
            sample_id = sample_id,
            bam_file = bam_file
    }

    output {
        File bai_file = SambambaIndex.bai_file
    }
}
