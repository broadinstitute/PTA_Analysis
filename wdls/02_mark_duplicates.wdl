version 1.0

# Workflow: MarkDuplicatesWorkflow
# Author: Shadi Zaheri
# Date: 2025-01-22
# Description: Marks duplicates in BAM files using GATK MarkDuplicates and creates BAM index using a separate task.

workflow MarkDuplicatesWorkflow {
  parameter_meta {
    input_bam: "Path to the input BAM file. Must be coordinate-sorted."
    sample_name: "Sample name used for output file naming."
    metrics_filename: "Path to the output duplication metrics file."
    preemptible_tries: "Number of preemptible retries allowed for the task (Google Cloud only)."
    memory_gb_markduplicates: "Memory allocated for the MarkDuplicates task in gigabytes."
    disk_gb_markduplicates: "Disk space allocated for the MarkDuplicates task in gigabytes."
    memory_gb_index: "Memory allocated for the IndexBam task in gigabytes."
    disk_gb_index: "Disk space allocated for the IndexBam task in gigabytes."
    cpu_markduplicates: "Number of CPU cores allocated for MarkDuplicates."
    cpu_index: "Number of CPU cores allocated for IndexBam."
    docker: "Docker image used to run the MarkDuplicates task."
  }

  input {
    File input_bam
    String sample_name
    String metrics_filename
    Int preemptible_tries = 1
    Int memory_gb_markduplicates = 16
    Int disk_gb_markduplicates = 200
    Int memory_gb_index = 4
    Int disk_gb_index = 50
    Int cpu_markduplicates = 4
    Int cpu_index = 1
    String docker = "broadinstitute/gatk:4.6.1.0"
  }

  # MarkDuplicates Step
  call MarkDuplicates {
    input:
      input_bam = input_bam,
      sample_name = sample_name,
      metrics_filename = metrics_filename,
      preemptible_tries = preemptible_tries,
      memory_gb = memory_gb_markduplicates,
      disk_gb = disk_gb_markduplicates,
      cpu = cpu_markduplicates,
      docker = docker
  }

  # IndexBam Step
  call IndexBam {
    input:
      input_bam = MarkDuplicates.output_bam,
      sample_name = sample_name,
      memory_gb = memory_gb_index,
      disk_gb = disk_gb_index,
      cpu = cpu_index,
      docker = docker
  }

  output {
    File output_bam = MarkDuplicates.output_bam
    File output_bam_index = IndexBam.output_bam_index
    File duplicate_metrics = MarkDuplicates.duplicate_metrics
  }
}

# MarkDuplicates Task
task MarkDuplicates {
  input {
    File input_bam
    String sample_name
    String metrics_filename
    Int preemptible_tries
    Int memory_gb
    Int disk_gb
    Int cpu
    String docker
  }

  command {
    set -euxo pipefail

    gatk MarkDuplicates \
      -I ~{input_bam} \
      -O ~{sample_name}.bam \
      -M ~{metrics_filename} \
      --TAGGING_POLICY All \
      --REMOVE_DUPLICATES false \
      --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
      --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 8000 \
      --MAX_OPTICAL_DUPLICATE_SET_SIZE 300000 \
      --MAX_RECORDS_IN_RAM 500000

    echo "MarkDuplicates completed successfully."
  }

  runtime {
    docker: docker
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_gb} HDD"
    cpu: cpu
    preemptible: preemptible_tries
  }

  output {
    File output_bam = "~{sample_name}.bam"
    File duplicate_metrics = "~{metrics_filename}"
  }
}

# IndexBam Task
task IndexBam {
  input {
    File input_bam
    String sample_name
    Int memory_gb
    Int disk_gb
    Int cpu
    String docker
  }

  command {
    set -euxo pipefail

    gatk BuildBamIndex \
      -I ~{input_bam} \
      -O ~{sample_name}.bam.bai

    echo "Indexing completed successfully."
  }

  runtime {
    docker: docker
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_gb} HDD"
    cpu: cpu
    preemptible: 1
  }

  output {
    File output_bam_index = "~{sample_name}.bam.bai"
  }
}
