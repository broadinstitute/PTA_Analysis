version 1.0

task tabix {
    input {
        String donor_id
        String sample_id
        File vcf_gz
        String? optional = "" # Making this optional parameter configurable
        
        # Runtime parameters
        String docker = "quay.io/biocontainers/htslib:1.15--h9753748_0"
        Int memory_gb = 4
        Int cpu = 1
        Int disk_size_gb = 50
    }

    command <<<
        set -euo pipefail
        hostname
        
        tabix ~{optional} ~{vcf_gz}
    >>>

    output {
        File tabix_index = "~{vcf_gz}.tbi"
    }

    runtime {
        docker: docker
        memory: "~{memory_gb}G"
        cpu: cpu
        disks: "local-disk ~{disk_size_gb} SSD"
        preemptible: 3
    }

    meta {
        author: "Converted from Nextflow"
        version: "1.0"
    }
}

workflow run_tabix {
    input {
        String donor_id
        String sample_id
        File vcf_gz
        String? optional
    }

    call tabix {
        input:
            donor_id = donor_id,
            sample_id = sample_id,
            vcf_gz = vcf_gz,
            optional = optional
    }

    output {
        File tabix_index = tabix.tabix_index
    }
}