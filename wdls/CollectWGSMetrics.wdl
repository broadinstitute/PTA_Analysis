version 1.0

# Workflow to run GATK CollectWGSMetrics using Terra inputs and customizable runtime parameters.
workflow CollectWGSMetricsWorkflow {

  parameter_meta {
    sample_id: "Unique sample identifier used to name the output metrics file."
    bam: "BAM file corresponding to the sample."
    ref_fasta: "Reference genome FASTA file."
    ref_index: "Reference genome FASTA index file."
    ref_dict: "Reference genome sequence dictionary file."
    optional: "Optional parameters to pass to the GATK CollectWGSMetrics command."
    cpu: "Number of CPUs allocated to the task."
    mem_gb: "Amount of memory in GB allocated to the task; Java heap is set to mem_gb - 4."
    disk_gb: "Amount of disk space (in GB) to allocate. This integer is converted into the runtime string 'local-disk <disk_gb> HDD'."
    preemptible: "Indicator for running the task on a preemptible node (0 for false, >0 for true)."
  }

  input {
    String sample_id
    File bam

    File ref_fasta
    File ref_index
    File ref_dict

    String optional = ""
    Int cpu = 1
    Int mem_gb = 8
    Int disk_gb = 200
    Int preemptible = 0
  }

  call CollectWGSMetrics {
    input:
      sample_id   = sample_id,
      bam         = bam,
      ref_fasta   = ref_fasta,
      ref_index   = ref_index,
      ref_dict    = ref_dict,
      optional    = optional,
      cpu         = cpu,
      mem_gb      = mem_gb,
      disk_gb     = disk_gb,
      preemptible = preemptible
  }

  output {
    File wgs_metrics = CollectWGSMetrics.wgs_metrics
  }
}

task CollectWGSMetrics {
  input {
    # Sample input: a sample id and its BAM file.
    String sample_id
    File bam

    # Reference inputs for Terra.
    File ref_fasta
    File ref_index
    File ref_dict

    # Additional command-line options (if any) to pass to GATK.
    String optional = ""

    # Runtime resource parameters.
    Int cpu = 1
    # "mem_gb" is used both for the runtime memory (in GB) and for computing Java's -Xmx value.
    Int mem_gb = 8
    # Disk space in gigabytes, which we then convert to the required string format.
    Int disk_gb = 180
    Int preemptible = 1
  }

  command <<<!
    gatk --java-options "-Xmx~{mem_gb - 4}g -Djava.io.tmpdir=$TMPDIR" CollectWgsMetrics \
      -I ~{bam} \
      -O ~{sample_id}.wgs_metrics.txt \
      -R ~{ref_fasta} \
      ~{optional}
    
    sed -i 's/picard\.analysis\.WgsMetrics/picard\.analysis\.CollectWGSMetrics$WgsMetrics/' ~{sample_id}.wgs_metrics.txt
  >>>

  output {
    File wgs_metrics = "~{sample_id}.wgs_metrics.txt"
  }

  runtime {
    docker: "broadinstitute/gatk:4.6.1.0"
    cpu: cpu
    memory: "~{mem_gb}GB"
    disk: "local-disk ~{disk_gb} HDD"
    preemptible: preemptible
  }
}