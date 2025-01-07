version 1.0

task GripssApplicationKt {
    input {
        String normal_sample_id
        String tumor_sample_id
        String donor_id
        File gridss_unfiltered_vcf
        File gridss_unfiltered_tbi
        File genome_fasta
        File breakpoint_hotspot
        File breakend_pod
        File breakpoint_pon
        String? additional_args
        
        # Runtime parameters
        Int memory_gb = 32
        Int cpu = 4
        String docker = "gridss/gridss-purple-linx:1.3.2"
    }

    Int java_mem_gb = memory_gb - 8

    command <<<
        java -Xmx~{java_mem_gb}g -cp /opt/hmftools/gripss-1.9.jar com.hartwig.hmftools.gripss.GripssApplicationKt \
        -ref_genome ~{genome_fasta} \
        -input_vcf ~{gridss_unfiltered_vcf} \
        -output_vcf ~{tumor_sample_id}.gripss.somatic.vcf.gz \
        -tumor ~{tumor_sample_id} \
        -reference ~{normal_sample_id} \
        -breakpoint_hotspot ~{breakpoint_hotspot} \
        -breakend_pon ~{breakend_pod} \
        -breakpoint_pon ~{breakpoint_pon} \
        ~{additional_args}
    >>>

    output {
        File gripss_somatic_vcf = "~{tumor_sample_id}.gripss.somatic.vcf.gz"
        File gripss_somatic_vcf_tbi = "~{tumor_sample_id}.gripss.somatic.vcf.gz.tbi"
    }

    runtime {
        docker: docker
        memory: "~{memory_gb} GB"
        cpu: cpu
    }
}

workflow Gripss {
    input {
        String normal_sample_id
        String tumor_sample_id
        String donor_id
        File gridss_unfiltered_vcf
        File gridss_unfiltered_tbi
        File genome_fasta
        File breakpoint_hotspot
        File breakend_pod
        File breakpoint_pon
        String? additional_args
    }

    call GripssApplicationKt {
        input:
            normal_sample_id = normal_sample_id,
            tumor_sample_id = tumor_sample_id,
            donor_id = donor_id,
            gridss_unfiltered_vcf = gridss_unfiltered_vcf,
            gridss_unfiltered_tbi = gridss_unfiltered_tbi,
            genome_fasta = genome_fasta,
            breakpoint_hotspot = breakpoint_hotspot,
            breakend_pod = breakend_pod,
            breakpoint_pon = breakpoint_pon,
            additional_args = additional_args
    }

    output {
        File somatic_vcf = GripssApplicationKt.gripss_somatic_vcf
        File somatic_vcf_tbi = GripssApplicationKt.gripss_somatic_vcf_tbi
    }
}