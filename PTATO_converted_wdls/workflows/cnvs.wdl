version 1.0

workflow cnvs {
  input {
    Array[File] bams
    Array[File] germline_vcfs
    Array[String] bulk_names
    String out_dir
    Boolean use_precomputed_cobalt = false
    String? cobalt_ratio_tsv_dir
  }

  scatter (i in range(length(bams))) {
    call classify_bams {
      input:
        bam_file = bams[i]
        bulk_name = bulk_names[i]
    }
  }

  # Step 1: Extract or generate cobalt ratio TSV files
  if (use_precomputed_cobalt) {
    call extractCobaltRatioTsvFromDir {
      input:
        cobalt_ratio_tsv_dir = cobalt_ratio_tsv_dir
    }
  } else {
    call get_cobalt_files {
      input:
        normal_bams = classify_bams.normal_bams
        tumor_bams = classify_bams.tumor_bams
    }
  }

  # Step 2: Filter CNV files
  call filter_cnv_files {
    input:
      cobalt_ratio_tsv_files = select_first([extractCobaltRatioTsvFromDir.cobalt_ratio_tsv_files, get_cobalt_files.cobalt_ratio_tsv_files])
      germline_vcfs = germline_vcfs
      normal_bams = classify_bams.normal_bams
      tumor_bams = classify_bams.tumor_bams
  }

  output {
    Array[File] filtered_cnv_files = filter_cnv_files.filtered_cnv_files
  }
}
task classify_bams {
  input {
    File bam_file
    String bulk_name
  }

  command <<<
    # Classify BAMs as normal or tumor based on bulk name
    if [[ "${bulk_name}" == "Normal" ]]; then
      echo "Normal BAM: ${bam_file}"
    else
      echo "Tumor BAM: ${bam_file}"
    fi
  >>>

  output {
    File classified_bam = bam_file
    String label = if "${bulk_name}" == "Normal" then "Normal" else "Tumor"
  }

  runtime {
    docker: "ubuntu:20.04"
    memory: "4 GB"
    cpu: 2
  }
}
task extractCobaltRatioTsvFromDir {
  input {
    String cobalt_ratio_tsv_dir
  }

  command <<<
    echo "Extracting cobalt ratio TSV files from ${cobalt_ratio_tsv_dir}"
    find ${cobalt_ratio_tsv_dir} -name "*.cobalt.ratio.tsv" > cobalt_files.list
  >>>

  output {
    Array[File] cobalt_ratio_tsv_files = read_lines("cobalt_files.list")
  }

  runtime {
    docker: "ubuntu:20.04"
    memory: "4 GB"
    cpu: 2
  }
}
task filter_cnv_files {
  input {
    Array[File] cobalt_ratio_tsv_files
    Array[File] germline_vcfs
    Array[File] normal_bams
    Array[File] tumor_bams
  }

  command <<<
    echo "Filtering CNV files..."
    # Mock command simulating CNV filtering logic
    for cobalt_file in ${sep=' ' cobalt_ratio_tsv_files}; do
      echo "Processing ${cobalt_file}"
    done
  >>>

  output {
    Array[File] filtered_cnv_files = glob("*.filtered.tsv")
  }

  runtime {
    docker: "cnv-filter:latest"
    memory: "8 GB"
    cpu: 4
  }
}
