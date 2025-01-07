version 1.0

task CountBamLinesApplication {
    input {
        String donor_id
        String normal_sample_id
        File normal_bam
        File normal_bai
        String tumor_sample_id
        File tumor_bam
        File tumor_bai
        File genome_fasta
        File gc_profile
        String? optional_args
        
        # Runtime parameters
        Int cpus = 4
        Int memory_gb = 16
        String docker = "gridss/gridss-purple-linx:1.3.2"
    }

    Int java_mem_gb = memory_gb - 4

    command <<<
        java -Xmx~{java_mem_gb}g \
        -Dsamjdk.reference_fasta=~{genome_fasta} \
        -Dsamjdk.use_async_io_read_samtools=true \
        -Dsamjdk.use_async_io_write_samtools=true \
        -Dsamjdk.use_async_io_write_tribble=true \
        -Dsamjdk.buffer_size=4194304 \
        -Dsamjdk.async_io_read_threads=~{cpus} \
        -cp /opt/hmftools/cobalt-1.11.jar com.hartwig.hmftools.cobalt.CountBamLinesApplication \
        -threads ~{cpus} \
        -tumor ~{tumor_sample_id} \
        -tumor_bam ~{tumor_bam} \
        -ref_genome ~{genome_fasta} \
        -output_dir ./~{tumor_sample_id} \
        -gc_profile ~{gc_profile} \
        -reference ~{normal_sample_id} \
        -reference_bam ~{normal_bam} \
        ~{optional_args}
    >>>

    output {
        Array[File] cobalt_files = glob("${tumor_sample_id}/*")
        File cobalt_ratio = "${tumor_sample_id}/${tumor_sample_id}.cobalt.ratio.tsv"
    }

    runtime {
        docker: docker
        cpu: cpus
        memory: "~{memory_gb}GB"
    }
}

workflow CountBamLinesWorkflow {
    input {
        String donor_id
        String normal_sample_id
        File normal_bam
        File normal_bai
        String tumor_sample_id
        File tumor_bam
        File tumor_bai
        File genome_fasta
        File gc_profile
        String? optional_args
    }

    call CountBamLinesApplication {
        input:
            donor_id = donor_id,
            normal_sample_id = normal_sample_id,
            normal_bam = normal_bam,
            normal_bai = normal_bai,
            tumor_sample_id = tumor_sample_id,
            tumor_bam = tumor_bam,
            tumor_bai = tumor_bai,
            genome_fasta = genome_fasta,
            gc_profile = gc_profile,
            optional_args = optional_args
    }

    output {
        Array[File] cobalt_files = CountBamLinesApplication.cobalt_files
        File cobalt_ratio = CountBamLinesApplication.cobalt_ratio
    }
}