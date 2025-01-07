version 1.0

workflow Shapeit {
    input {
        String donor_id
        String sample_id
        File vcf_gz
        File tbi
        String chrom
        File reference_vcf
        File genetic_map
        Int threads = 8
        String? optional_args = ""
    }

    call ShapeitPhasing {
        input:
            donor_id = donor_id,
            sample_id = sample_id,
            vcf_gz = vcf_gz,
            tbi = tbi,
            chrom = chrom,
            reference_vcf = reference_vcf,
            genetic_map = genetic_map,
            threads = threads,
            optional_args = optional_args
    }

    output {
        File phased_vcf = ShapeitPhasing.phased_vcf
    }
}

task ShapeitPhasing {
    input {
        String donor_id
        String sample_id
        File vcf_gz
        File tbi
        String chrom
        File reference_vcf
        File genetic_map
        Int threads
        String? optional_args
    }

    String output_filename = "~{sample_id}.~{chrom}.phased.vcf.gz"

    command <<<
        set -euo pipefail
        hostname

        shapeit4 \
        --input ~{vcf_gz} \
        --region ~{chrom} \
        --output ~{output_filename} \
        --map ~{genetic_map} \
        --thread ~{threads} \
        --sequencing \
        --reference ~{reference_vcf} \
        ~{optional_args}
    >>>

    runtime {
        docker: "quay.io/biocontainers/shapeit4:4.2.2--h24bf969_1"
        cpu: threads
        memory: "16 GB"
    }

    output {
        File phased_vcf = output_filename
    }
}