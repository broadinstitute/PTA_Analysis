version 1.0

workflow MarkDuplicatesWorkflow {
  input {
    File input_bam
    String output_bam_basename
    String metrics_filename
    Int preemptible_tries = 1
    Int memory_gb = 8
    Int disk_gb = 50
    Int cpu = 4
    String docker = "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
  }

  call MarkDuplicates {
    input:
      input_bam = input_bam,
      output_bam_basename = output_bam_basename,
      metrics_filename = metrics_filename,
      preemptible_tries = preemptible_tries,
      memory_gb = memory_gb,
      disk_gb = disk_gb,
      cpu = cpu,
      docker = docker
  }

  output {
    File output_bam = MarkDuplicates.output_bam
    File output_bam_index = MarkDuplicates.output_bam_index
    File duplicate_metrics = MarkDuplicates.duplicate_metrics
  }
}

task MarkDuplicates {
  input {
    File input_bam
    String output_bam_basename
    String metrics_filename
    Int preemptible_tries
    Int memory_gb
    Int disk_gb
    Int cpu
    String docker
  }

  command {
    java -Xmx~{memory_gb}G -jar /usr/picard/picard.jar MarkDuplicates \
      INPUT=~{input_bam} \
      OUTPUT=~{output_bam_basename}.bam \
      METRICS_FILE=~{metrics_filename} \
      CREATE_INDEX=true
  }

  runtime {
    docker: docker
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_gb} HDD"
    cpu: cpu  # Static integer value
    preemptible: preemptible_tries
  }

  output {
    File output_bam = "~{output_bam_basename}.bam"
    File output_bam_index = "~{output_bam_basename}.bam.bai"
    File duplicate_metrics = "~{metrics_filename}"
  }
}
