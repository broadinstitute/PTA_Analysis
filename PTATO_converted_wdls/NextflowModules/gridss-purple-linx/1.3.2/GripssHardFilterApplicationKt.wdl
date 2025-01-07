version 1.0

task GripssHardFilterApplicationKt {
    input {
        String donor_id
        String normal_sample_id
        String tumor_sample_id
        File gripss_somatic_vcf
        File gripss_somatic_tbi
        
        # Runtime parameters
        String docker = "gridss/gridss-purple-linx:1.3.2"
        Int memory_gb = 16
        Int disk_size_gb = 50
    }

    command <<<
        java -Xmx~{memory_gb-4}g -cp /opt/hmftools/gripss-1.9.jar \
        com.hartwig.hmftools.gripss.GripssHardFilterApplicationKt \
        -input_vcf ~{gripss_somatic_vcf} \
        -output_vcf ~{tumor_sample_id}.gripss.somatic.filtered.vcf.gz
    >>>

    output {
        File filtered_vcf = "~{tumor_sample_id}.gripss.somatic.filtered.vcf.gz"
        File filtered_vcf_index = "~{tumor_sample_id}.gripss.somatic.filtered.vcf.gz.tbi"
    }

    runtime {
        docker: docker
        memory: "~{memory_gb}GB"
        disks: "local-disk ~{disk_size_gb} SSD"
        preemptible: 3
    }
}

workflow GripssHardFilter {
    input {
        String donor_id
        String normal_sample_id
        String tumor_sample_id
        File gripss_somatic_vcf
        File gripss_somatic_tbi
    }

    call GripssHardFilterApplicationKt {
        input:
            donor_id = donor_id,
            normal_sample_id = normal_sample_id,
            tumor_sample_id = tumor_sample_id,
            gripss_somatic_vcf = gripss_somatic_vcf,
            gripss_somatic_tbi = gripss_somatic_tbi
    }

    output {
        File filtered_vcf = GripssHardFilterApplicationKt.filtered_vcf
        File filtered_vcf_index = GripssHardFilterApplicationKt.filtered_vcf_index
    }
}