version 1.0

workflow filter_cnv_files {
  input {
    Array[File] cobalt_ratio_tsv_files
    Array[File] germline_vcf_files
    Array[File] normal_bams
    Array[File] tumor_bams
    Map[String, String] optional_cnvs_dirs
    String out_dir
  }

  # Step 1: Extract or Filter Cobalt files
  if (defined(optional_cnvs_dirs["cobalt_filtered_readcounts_dir"])) {
    call ExtractCobaltFilteredReadCounts {
      input:
        cobalt_filtered_readcounts_dir = optional_cnvs_dirs["cobalt_filtered_readcounts_dir"]
    }
    call ExtractCobaltFilteredReadCountsSegments {
      input:
        cobalt_filtered_readcounts_dir = optional_cnvs_dirs["cobalt_filtered_readcounts_dir"]
    }
  } else {
    call FilterCobalt {
      input:
        cobalt_ratio_tsv_files = cobalt_ratio_tsv_files
    }
  }

  # Step 2: Extract or Filter BAF files
  if (defined(optional_cnvs_dirs["baf_filtered_files_dir"])) {
    call ExtractBafFilteredFiles {
      input:
        baf_filtered_files_dir = optional_cnvs_dirs["baf_filtered_files_dir"]
    }
    call ExtractBafSegmentsFiles {
      input:
        baf_filtered_files_dir = optional_cnvs_dirs["baf_filtered_files_dir"]
    }
    call ExtractBafBinnedFiles {
      input:
        baf_filtered_files_dir = optional_cnvs_dirs["baf_filtered_files_dir"]
    }
  } else {
    scatter (sample_id in normal_bams) {
      call FilterBAF {
        input:
          normal_bam = sample_id
          tumor_bam = tumor_bam
      }
    }
  }

  output {
    Array[File] filtered_cnv_files
  }
}
