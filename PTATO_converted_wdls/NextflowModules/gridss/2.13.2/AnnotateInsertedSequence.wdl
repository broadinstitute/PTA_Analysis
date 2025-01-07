version 1.0

task AnnotateInsertedSequence {
    input {
        String donor_id
        String normal_sample_id
        String tumor_sample_id
        File gridss_driver_vcf
        File gridss_driver_tbi
        File viral_reference
        Int memory_gb = 16
        Int cpus = 8
        String docker = "gridss/gridss:2.13.2"
    }

    Int java_memory_gb = memory_gb - 4

    command <<<
        java -Xmx~{java_memory_gb}g \
        -Dsamjdk.create_index=true \
        -Dsamjdk.use_async_io_read_samtools=true \
        -Dsamjdk.use_async_io_write_samtools=true \
        -Dsamjdk.use_async_io_write_tribble=true \
        -Dsamjdk.buffer_size=4194304 \
        -cp /opt/gridss/gridss-2.13.2-gridss-jar-with-dependencies.jar gridss.AnnotateInsertedSequence \
        REFERENCE_SEQUENCE=~{viral_reference} \
        INPUT=~{gridss_driver_vcf} \
        OUTPUT=~{tumor_sample_id}.gridss.unfiltered.vcf.gz \
        ALIGNMENT=APPEND \
        WORKER_THREADS=~{cpus}
    >>>

    output {
        File gridss_unfiltered_vcf = "~{tumor_sample_id}.gridss.unfiltered.vcf.gz"
        File gridss_unfiltered_vcf_tbi = "~{tumor_sample_id}.gridss.unfiltered.vcf.gz.tbi"
    }

    runtime {
        docker: docker
        memory: "~{memory_gb}GB"
        cpu: cpus
    }
}

workflow AnnotateInsertedSequenceWorkflow {
    input {
        String donor_id
        String normal_sample_id
        String tumor_sample_id
        File gridss_driver_vcf
        File gridss_driver_tbi
        File viral_reference
    }

    call AnnotateInsertedSequence {
        input:
            donor_id = donor_id,
            normal_sample_id = normal_sample_id,
            tumor_sample_id = tumor_sample_id,
            gridss_driver_vcf = gridss_driver_vcf,
            gridss_driver_tbi = gridss_driver_tbi,
            viral_reference = viral_reference
    }

    output {
        File gridss_unfiltered_vcf = AnnotateInsertedSequence.gridss_unfiltered_vcf
        File gridss_unfiltered_vcf_tbi = AnnotateInsertedSequence.gridss_unfiltered_vcf_tbi
    }
}
