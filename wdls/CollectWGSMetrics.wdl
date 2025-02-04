version 1.0

#     description: "This task runs GATK's CollectWgsMetrics to compute whole-genome sequencing metrics."
#    documentation: "https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-CollectWgsMetrics"
workflow CollectWGSMetricsWorkflow {
  
  input {
    File input_bam
    File input_bai
    File reference_fasta
    File reference_fasta_index
    File reference_dict
    Int coverage_cap = 10000
    String output_prefix = "wgs_metrics"

    # Terra-configurable runtime options
    Int cpu = 2
    String memory = "8G"
    String disk_size = "50G"
    Int preemptible_tries = 3  # Default: 3 retries on preemptible instances
  }

  parameter_meta {
    input_bam: "The input BAM file containing aligned reads."
    input_bai: "The index file for the input BAM."
    reference_fasta: "Reference genome in FASTA format."
    reference_fasta_index: "Index file for the reference genome."
    reference_dict: "Sequence dictionary file for the reference genome."
    coverage_cap: "The maximum coverage depth considered in calculations to avoid extreme outliers. Default is 10,000."
    output_prefix: "Prefix for the output metrics file."
    cpu: "Number of CPU cores allocated for the task."
    memory: "Memory allocated for the task."
    disk_size: "Disk size allocated for the task."
    preemptible_tries: "Number of preemptible retries allowed for each task."
  }

  call CollectWGSMetricsTask {
    input:
      input_bam = input_bam,
      input_bai = input_bai,
      reference_fasta = reference_fasta,
      reference_fasta_index = reference_fasta_index,
      reference_dict = reference_dict,
      coverage_cap = coverage_cap,
      output_prefix = output_prefix,
      cpu = cpu,
      memory = memory,
      disk_size = disk_size,
      preemptible_tries = preemptible_tries
  }

  output {
    File wgs_metrics = CollectWGSMetricsTask.wgs_metrics
  }
}

task CollectWGSMetricsTask {
  input {
    File input_bam
    File input_bai
    File reference_fasta
    File reference_fasta_index
    File reference_dict
    Int coverage_cap
    String output_prefix

    # Terra-configurable runtime options
    Int cpu
    String memory
    String disk_size
    Int preemptible_tries
  }

  command {
    gatk CollectWgsMetrics \
      -I ~{input_bam} \
      -O ~{output_prefix}.txt \
      -R ~{reference_fasta} \
      --COVERAGE_CAP ~{coverage_cap}
  }

  output {
    File wgs_metrics = "~{output_prefix}.txt"
  }

  runtime {
    docker: "broadinstitute/gatk:4.6.1.0"#"broadinstitute/gatk:latest"
    cpu: cpu
    memory: memory
    disks: "local-disk ~{disk_size} HDD"
    preemptible_tries: preemptible_tries  # No conversion needed!
  }
}
