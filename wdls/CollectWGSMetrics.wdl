version 1.0

# Summary:
# This workflow runs GATK CollectWgsMetrics to compute whole-genome sequencing metrics on a given BAM file.
# It outputs a metrics file containing coverage statistics and other quality metrics.

workflow CollectWGSMetricsWorkflow {
    
    input {
        String sample_id
        File bam
        File genome_fasta
        String optional_parameters
        Int memory_gb
    }
    
    call CollectWGSMetrics {
        input:
            sample_id = sample_id,
            bam = bam,
            genome_fasta = genome_fasta,
            optional_parameters = optional_parameters,
            memory_gb = memory_gb
    }
    
    output {
        File wgs_metrics = CollectWGSMetrics.wgs_metrics
    }
}

task CollectWGSMetrics {
    
    input {
        String sample_id # Unique identifier for the sample
        File bam # Input BAM file for whole-genome sequencing analysis
        File genome_fasta # Reference genome FASTA file
        String optional_parameters # Additional optional parameters for GATK CollectWgsMetrics
        Int memory_gb # Memory allocation for the task (in GB)
    }
    
    meta {
        description: "Runs GATK CollectWgsMetrics to generate whole-genome sequencing quality metrics."
    }
    
    command {
        gatk --java-options "-Xmx~{memory_gb - 4}G -Djava.io.tmpdir=$TMPDIR" \
        CollectWgsMetrics \
        -I ~{bam} \
        -O ~{sample_id}.wgs_metrics.txt \
        -R ~{genome_fasta} \
        ~{optional_parameters}
        
        sed -i 's/picard\.analysis\.WgsMetrics/picard\.analysis\.CollectWgsMetrics\$WgsMetrics/' ~{sample_id}.wgs_metrics.txt
    }
    
    runtime {
        docker: "broadinstitute/gatk:4.6.1.0"
        memory: "~{memory_gb} GB"
        cpu: 1
    }
    
    output {
        File wgs_metrics = "~{sample_id}.wgs_metrics.txt"
    }
}
