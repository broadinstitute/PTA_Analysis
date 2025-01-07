version 1.0

task bgzip {
    input {
        String donor_id
        String sample_id
        File vcf
        String? optional_args = ""
        String docker = "quay.io/biocontainers/htslib:1.15--h9753748_0"
    }

    command <<<
        hostname
        bgzip -c ~{vcf} ~{optional_args} > ~{vcf}.gz
    >>>

    output {
        File vcf_gz = "~{vcf}.gz"
    }

    runtime {
        docker: docker
        memory: "4 GB"  # You can adjust these runtime parameters
        cpu: 1
        disks: "local-disk 50 SSD"
    }
}

workflow bgzip_workflow {
    input {
        String donor_id
        String sample_id
        File vcf
    }

    call bgzip {
        input:
            donor_id = donor_id,
            sample_id = sample_id,
            vcf = vcf
    }

    output {
        File output_vcf_gz = bgzip.vcf_gz
    }
}