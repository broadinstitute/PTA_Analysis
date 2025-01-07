version 1.0

workflow bedtools_sort_workflow {
    input {
        String donor_id
        String sample_id
        File bed_file
        String? optional_args
    }

    call bedtools_sort {
        input:
            donor_id = donor_id,
            sample_id = sample_id,
            bed_file = bed_file,
            optional_args = optional_args
    }

    output {
        File sorted_bed = bedtools_sort.sorted_bed
    }

    meta {
        author: "Converted from Nextflow PTATO by Shadi Zaheri"
        version: "1.0.0"
        description: "Workflow for sorting BED files using bedtools sort"
    }
}

task bedtools_sort {
    input {
        String donor_id
        String sample_id
        File bed_file
        String? optional_args = ""  # Equivalent to params.bedtoolssort.optional
        
        # Runtime parameters
        String docker = "quay.io/biocontainers/bedtools:2.30.0--h468198e_3"
        Int? memory_gb = 4
        Int? disk_size_gb = 50
    }

    command <<<
        set -euo pipefail
        hostname

        bedtools sort \
            -i ~{bed_file} \
            ~{optional_args} \
            > ~{sample_id}.sorted.bed
    >>>

    output {
        File sorted_bed = "~{sample_id}.sorted.bed"
    }

    runtime {
        docker: docker
        memory: "~{memory_gb}G"
        disks: "local-disk ~{disk_size_gb} SSD"
        maxRetries: 3
    }

    meta {
        author: "Converted from Nextflow PTATO by Shadi Zaheri"
        version: "2.30.0"
    }
}
