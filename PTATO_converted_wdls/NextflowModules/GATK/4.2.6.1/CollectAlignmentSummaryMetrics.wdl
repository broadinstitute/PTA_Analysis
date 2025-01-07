version 1.0

workflow CollectAlignmentSummaryMetricsWorkflow {
    input {
        String donor_id
        String sample_id
        File input_bam
        File input_bai
        File ref_fasta
        String? optional_args
        Int memory_gb = 16
        Int disk_size_gb = 100
    }

    call CollectAlignmentSummaryMetrics {
        input:
            donor_id = donor_id,
            sample_id = sample_id,
            input_bam = input_bam,
            input_bai = input_bai,
            ref_fasta = ref_fasta,
            optional_args = optional_args,
            memory_gb = memory_gb,
            disk_size_gb = disk_size_gb
    }

    output {
        File alignment_summary_metrics = CollectAlignmentSummaryMetrics.metrics_output
    }
}

task CollectAlignmentSummaryMetrics {
    input {
        String donor_id
        String sample_id
        File input_bam
        File input_bai
        File ref_fasta
        String? optional_args
        Int memory_gb
        Int disk_size_gb
    }

    Int java_memory_gb = memory_gb - 4

    command <<<
        set -euo pipefail

        gatk --java-options "-Xmx~{java_memory_gb}G" \
            CollectAlignmentSummaryMetrics \
            -I ~{input_bam} \
            -O ~{sample_id}.alignment_summary_metrics.txt \
            -R ~{ref_fasta} \
            ~{optional_args}

        sed -i 's/picard\.analysis\.AlignmentSummaryMetrics/picard\.analysis\.CollectAlignmentSummaryMetrics\$AlignmentSummaryMetrics/' ~{sample_id}.alignment_summary_metrics.txt
    >>>

    output {
        File metrics_output = "~{sample_id}.alignment_summary_metrics.txt"
    }

    runtime {
        docker: "broadinstitute/gatk:4.2.6.1"
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_size_gb} HDD"
        preemptible: 3
    }
}
