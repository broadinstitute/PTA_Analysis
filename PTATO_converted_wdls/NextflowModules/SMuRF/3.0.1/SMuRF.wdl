version 1.0

task smurf {
    input {
        String donor_id
        String germline_sample_id
        File germline_vcf
        File germline_tbi
        Array[String] bam_sample_ids
        Array[File] bam_files
        Array[File] bai_files
        Array[String] bulk_names
        String smurf_config
        Int cpus = 4
        String docker = "vanboxtelbioinformatics/smurf:3.0.1"
    }

    command <<<
        set -euo pipefail
        
        hostname

        python /smurf/SMuRF.py \
        -i ~{germline_vcf} \
        ~{if length(bam_files) > 0 then " -b " + sep(" -b ", bam_files) else ""} \
        ~{if length(bulk_names) > 0 then " -n " + sep(" -n ", bulk_names) else ""} \
        -t ~{cpus} \
        -c ~{smurf_config}

        bash /smurf/scripts/split_in_single_sample_vcfs.sh ~{germline_sample_id}.SMuRF.filtered.vcf

        # Remove specific bulk VCF files
        for BULK in ~{sep(" ", bulk_names)}; do
            if [[ "${BULK}" != "-n" ]]; then
                rm ~{germline_sample_id}_${BULK}.SMuRF.filtered.vcf*
            fi
        done
    >>>

    output {
        File smurf_vcf = "~{germline_sample_id}.SMuRF.vcf"
        File filtered_vcf = "~{germline_sample_id}.SMuRF.filtered.vcf"
        File vaf_plot = "~{germline_sample_id}.SMuRF.vafplot.pdf"
        Array[File] single_sample_vcfs = glob("~{germline_sample_id}_*.SMuRF.filtered.vcf")
    }

    runtime {
        docker: docker
        cpu: cpus
        memory: "16 GB"  # You may want to adjust this based on your needs
    }

    meta {
        author: "Converted from Nextflow"
        description: "SMuRF somatic variant calling task"
    }
}

workflow SMuRF {
    input {
        String donor_id
        String germline_sample_id
        File germline_vcf
        File germline_tbi
        Array[String] bam_sample_ids
        Array[File] bam_files
        Array[File] bai_files
        Array[String] bulk_names
        String smurf_config
    }

    call smurf {
        input:
            donor_id = donor_id,
            germline_sample_id = germline_sample_id,
            germline_vcf = germline_vcf,
            germline_tbi = germline_tbi,
            bam_sample_ids = bam_sample_ids,
            bam_files = bam_files,
            bai_files = bai_files,
            bulk_names = bulk_names,
            smurf_config = smurf_config
    }

    output {
        File smurf_vcf = smurf.smurf_vcf
        File filtered_vcf = smurf.filtered_vcf
        File vaf_plot = smurf.vaf_plot
        Array[File] single_sample_vcfs = smurf.single_sample_vcfs
    }
}