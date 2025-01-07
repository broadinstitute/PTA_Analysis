version 1.0

task gridss {
    input {
        String donor_id
        String normal_sample_id
        File normal_bam
        File normal_bai
        String tumor_sample_id
        File tumor_bam
        File tumor_bai
        File genome_fasta
        String? gridss_optional
        Int cpu = 8
        Int memory_gb = 32
    }

    command <<<
        gridss \
        --jvmheap ~{memory_gb-4}g \
        -o ~{tumor_sample_id}.gridss.driver.vcf.gz \
        -r ~{genome_fasta} \
        -t ~{cpu} \
        --labels ~{normal_sample_id},~{tumor_sample_id} \
        ~{default="" gridss_optional} \
        ~{normal_bam} \
        ~{tumor_bam}
    >>>

    output {
        File gridss_vcf = "~{tumor_sample_id}.gridss.driver.vcf.gz"
        File gridss_vcf_index = "~{tumor_sample_id}.gridss.driver.vcf.gz.tbi"
        File gridss_assembly_bam = "~{tumor_sample_id}.gridss.driver.vcf.gz.assembly.bam"
    }

    runtime {
        docker: "gridss/gridss:2.13.2"
        cpu: cpu
        memory: "~{memory_gb}GB"
    }
}

workflow run_gridss {
    input {
        String donor_id
        String normal_sample_id
        File normal_bam
        File normal_bai
        String tumor_sample_id
        File tumor_bam
        File tumor_bai
        File genome_fasta
        String? gridss_optional
    }

    call gridss {
        input:
            donor_id = donor_id,
            normal_sample_id = normal_sample_id,
            normal_bam = normal_bam,
            normal_bai = normal_bai,
            tumor_sample_id = tumor_sample_id,
            tumor_bam = tumor_bam,
            tumor_bai = tumor_bai,
            genome_fasta = genome_fasta,
            gridss_optional = gridss_optional
    }

    output {
        File gridss_vcf = gridss.gridss_vcf
        File gridss_vcf_index = gridss.gridss_vcf_index
        File gridss_assembly_bam = gridss.gridss_assembly_bam
    }
}
