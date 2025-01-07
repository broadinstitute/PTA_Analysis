version 1.0

workflow SMuRF_workflow {
    input {
        String donor_id
        String germline_sample_id
        File germline_vcf
        File germline_tbi
        Array[String] bam_sample_ids
        Array[File] bam_files
        Array[File] bai_files
        Array[String] bulk_names
        File smurf_config
        Int cpus = 4
        String docker = "vanboxtelbioinformatics/smurf:3.0.2"
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
            smurf_config = smurf_config,
            cpus = cpus,
            docker = docker
    }

    output {
        File smurf_vcf = smurf.smurf_vcf
        File filtered_vcf = smurf.filtered_vcf
        File vaf_plot = smurf.vaf_plot
        Array[File] single_sample_vcfs = smurf.single_sample_vcfs
    }
}

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
        File smurf_config
        String docker
        Int cpus
    }

    command <<<
        set -euo pipefail
        
        hostname

        python /smurf/SMuRF.py \
        -i ~{germline_vcf} \
        ~{sep=" -b " bam_files} \
        ~{sep=" -n " bulk_names} \
        -t ~{cpus} \
        -c ~{smurf_config}

        bash /smurf/scripts/split_in_single_sample_vcfs.sh ~{germline_sample_id}.SMuRF.filtered.vcf

        # Remove specific bulk VCFs
        for BULK in ~{sep=" " bulk_names}; do
            rm ~{germline_sample_id}_${BULK}.SMuRF.filtered.vcf*
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
        memory: "16 GB"
        disks: "local-disk 100 SSD"
        maxRetries: 3
    }
}