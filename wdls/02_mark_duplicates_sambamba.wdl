version 1.0

workflow MarkDuplicatesWithSambamba {
  input {
    File input_bam
    String output_bam_basename
    Int preemptible_tries = 1
    Int memory_gb = 8
    Int disk_gb = 50
    String sambamba_docker = "biocontainers/sambamba:v0.6.8_cv1"
  }

  call SambambaMarkDuplicates {
    input:
      input_bam = input_bam,
      output_bam_basename = output_bam_basename,
      preemptible_tries = preemptible_tries,
      memory_gb = memory_gb,
      disk_gb = disk_gb,
      sambamba_docker = sambamba_docker
  }

  output {
    File output_bam = SambambaMarkDuplicates.output_bam
    File output_bam_index = SambambaMarkDuplicates.output_bam_index
  }
}

task SambambaMarkDuplicates {
  input {
    File input_bam
    String output_bam_basename
    Int preemptible_tries
    Int memory_gb
    Int disk_gb
    String sambamba_docker
  }

  command {
    sambamba markdup \
      --nthreads=4 \
      --tmpdir=/tmp \
      --overflow-list-size 1000000 \
      ~{input_bam} ~{output_bam_basename}.bam
    sambamba index ~{output_bam_basename}.bam
  }

  runtime {
    docker: sambamba_docker
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: preemptible_tries
    cpu: 4
  }

  output {
    File output_bam = "~{output_bam_basename}.bam"
    File output_bam_index = "~{output_bam_basename}.bam.bai"
  }
}