version 2.2.0

workflow WalkerWorkflow {
  input {
    String donor_id
    String germline_sample_id
    File germline_vcf
    File germline_tbi
    Array[String] bam_sample_ids
    Array[File] bam_files
    Array[File] bai_files
    String sample_id
    File somatic_vcf
    File somatic_tbi
    Int cpu = 4
  }

  call Walker {
    input:
      donor_id = donor_id,
      germline_sample_id = germline_sample_id,
      germline_vcf = germline_vcf,
      germline_tbi = germline_tbi,
      bam_sample_ids = bam_sample_ids,
      bam_files = bam_files,
      bai_files = bai_files,
      sample_id = sample_id,
      somatic_vcf = somatic_vcf,
      somatic_tbi = somatic_tbi,
      cpu = cpu
  }

  output {
    File walker_vcf = Walker.walker_vcf
    File walker_bed = Walker.walker_bed
    File walker_txt = Walker.walker_txt
  }
}


task Walker {
  input {
    String donor_id
    String germline_sample_id
    File germline_vcf
    File germline_tbi
    Array[String] bam_sample_ids
    Array[File] bam_files
    Array[File] bai_files
    String sample_id
    File somatic_vcf
    File somatic_tbi
    Int cpu
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    b=""
    if [[ ! -z "~{sep=' ' bam_files}" ]]; then
      b=$(printf " -b %s" "~{sep=' -b ' bam_files}")
    fi

    python /walker/walker.py \
      -g ~{germline_vcf} \
      -s ~{somatic_vcf} \
      -t ~{cpu} \
      ${b} \
      -o ~{sample_id} \
      -f vcf -f bed -f txt
  >>>

  output {
    File walker_vcf = "~{sample_id}.walker.vcf"
    File walker_bed = "~{sample_id}.walker.bed"
    File walker_txt = "~{sample_id}.walker.txt"
  }

  runtime {
    docker: "vanboxtelbioinformatics/walker:2.2.0"
    cpu: ~{cpu}
    memory: "8 GB"
  }
}
