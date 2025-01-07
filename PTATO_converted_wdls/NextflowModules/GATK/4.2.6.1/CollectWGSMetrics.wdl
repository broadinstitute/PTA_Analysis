version 1.0

workflow CollectWGSMetricsWorkflow {
    input {
        String donor_id
        String sample_id
        File input_bam
        File input_bai
        File ref_fasta
        String? optional_args
        Int memory_gb = 16
        String docker = "broadinstitute/gatk:4.2.6.1"
    }

    call CollectWGSMetrics {
        input:
            donor_id = donor_id,
            sample_id = sample_id,
            input_bam = input_bam,
            input_bai = input_bai,
            ref_fasta = ref_fasta,
            optional_args = optional_args,
            memory_gb = memory_gb,
            docker = docker
    }

    output {
        File wgs_metrics = CollectWGSMetrics.output_metrics
    }
}

task CollectWGSMetrics {
    input {
        String donor_id
        String sample_id
        File input_bam
        File input_bai
        File ref_fasta
        String? optional_args
        Int memory_gb
        String docker
    }

    Int command_mem_gb = memory_gb - 4

    command <<<
        set -euo pipefail
        
        gatk --java-options "-Xmx~{command_mem_gb}G" \
        CollectWgsMetrics \
        -I ~{input_bam} \
        -O ~{sample_id}.wgs_metrics.txt \
        -R ~{ref_fasta} \
        ~{optional_args}

        sed -i 's/picard\.analysis\.WgsMetrics/picard\.analysis\.CollectWgsMetrics\$WgsMetrics/' ~{sample_id}.wgs_metrics.txt
    >>>

    output {
        File output_metrics = "~{sample_id}.wgs_metrics.txt"
    }

    runtime {
        docker: docker
        memory: "~{memory_gb} GB"
        disks: "local-disk 50 SSD"
    }
}
