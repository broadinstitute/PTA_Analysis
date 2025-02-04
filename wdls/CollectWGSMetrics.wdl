version 1.0

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
    # "mem_gb" is used both for the runtime memory (in GB) and for computing Java’s -Xmx value.
    Int mem_gb = 8
    String disk = "local-disk 50 HDD"
    Int preemptible = 0
  }

  command <<<
    gatk --java-options "-Xmx~{mem_gb - 4}g -Djava.io.tmpdir=$TMPDIR" CollectWgsMetrics \
      -I ~{bam} \
      -O ~{sample_id}.wgs_metrics.txt \
      -R ~{ref_fasta} \
      ~{optional}
    
    sed -i 's/picard\\.analysis\\.WgsMetrics/picard\\.analysis\\.CollectWGSMetrics\\\$WgsMetrics/' ~{sample_id}.wgs_metrics.txt
  >>>

  output {
    File wgs_metrics = "~{sample_id}.wgs_metrics.txt"
  }

  runtime {
    docker: "library://sawibo/default/bioinf-tools:gatk4.1.3.0"
    cpu: cpu
    memory: "~{mem_gb}GB"
    disk: disk
    preemptible: preemptible
  }

  meta {
    description: "Task to run GATK CollectWGSMetrics on a BAM file using Terra reference inputs."
  }

  parameter_meta {
    sample_id: {
      description: "Unique identifier for the sample, used to name the output file."
    }
    bam: {
      description: "Input BAM file for the sample."
    }
    ref_fasta: {
      description: "Reference genome FASTA file."
    }
    ref_index: {
      description: "Index file for the reference genome FASTA."
    }
    ref_dict: {
      description: "Sequence dictionary for the reference genome."
    }
    optional: {
      description: "Optional additional parameters to pass to the GATK command."
    }
    cpu: {
      description: "Number of CPUs to allocate for the task."
    }
    mem_gb: {
      description: "Total memory in GB to allocate for the task; Java heap is set to mem_gb - 4."
    }
    disk: {
      description: "Disk resource specification."
    }
    preemptible: {
      description: "Preemptible flag for running the task on preemptible nodes (0 = false, >0 = true)."
    }
  }
}

workflow CollectWGSMetricsWorkflow {
  input {
    String sample_id
    File bam

    File ref_fasta
    File ref_index
    File ref_dict

    String optional = ""
    Int cpu = 1
    Int mem_gb = 8
    String disk = "local-disk 50 HDD"
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
      disk        = disk,
      preemptible = preemptible
  }

  output {
    File wgs_metrics = CollectWGSMetrics.wgs_metrics
  }

  meta {
    description: "Workflow to run GATK CollectWGSMetrics using Terra inputs and customizable runtime parameters."
  }

  parameter_meta {
    sample_id: {
      description: "Unique sample identifier used to name the output metrics file."
    }
    bam: {
      description: "BAM file corresponding to the sample."
    }
    ref_fasta: {
      description: "Reference genome FASTA file."
    }
    ref_index: {
      description: "Reference genome FASTA index file."
    }
    ref_dict: {
      description: "Reference genome sequence dictionary file."
    }
    optional: {
      description: "Optional parameters to pass to the GATK CollectWGSMetrics command."
    }
    cpu: {
      description: "Number of CPUs allocated to the task."
    }
    mem_gb: {
      description: "Amount of memory in GB allocated to the task; Java heap is set to mem_gb - 4."
    }
    disk: {
      description: "Disk specification for the runtime."
    }
    preemptible: {
      description: "Indicator for running the task on a preemptible node (0 for false, >0 for true)."
    }
  }
}
