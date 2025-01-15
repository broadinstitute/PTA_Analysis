version 1.0

workflow MarkDuplicatesWorkflow {
  input {
    File input_bam
    String output_bam_basename
    String metrics_filename
    Int preemptible_tries = 1
    String docker = "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
  }

  call MarkDuplicates {
    input:
      input_bam = input_bam,
      output_bam_basename = output_bam_basename,
      metrics_filename = metrics_filename,
      preemptible_tries = preemptible_tries,
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
    String docker
  }

  command {
    java -Xmx4G -jar /usr/picard/picard.jar MarkDuplicates \
      INPUT=~{input_bam} \
      OUTPUT=~{output_bam_basename}.bam \
      METRICS_FILE=~{metrics_filename} \
      CREATE_INDEX=true
  }

  runtime {
    docker: docker
    memory: "8 GiB"
    disks: "local-disk 50 HDD"
    preemptible: preemptible_tries
  }

  output {
    File output_bam = "~{output_bam_basename}.bam"
    File output_bam_index = "~{output_bam_basename}.bam.bai"
    File duplicate_metrics = "~{metrics_filename}"
  }
}
