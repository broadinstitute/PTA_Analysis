version 1.0

workflow ptato {
  input {
    Array[String] bulk_names
    String input_vcfs_dir
    String bams_dir
    Boolean run_QC
    Boolean run_postqc
    Boolean run_snvs
    Boolean run_indels
    Boolean run_svs
    Boolean run_cnvs
    Map[String, String] optional_postqc_dirs
  }

  # Step 1: Generate unique donor IDs
  scatter (bulk_name in bulk_names) {
    call GenerateDonorIds {
      input:
        bulk_name = bulk_name
    }
  }

  # Step 2: Extract input VCFs and BAMs
  call ExtractInputVcfFromDir {
    input:
      input_vcfs_dir = input_vcfs_dir
  }

  call ExtractBamsFromDir {
    input:
      bams_dir = bams_dir
  }

  # Step 3: Process BAM files (indexing)
  call GetIndexedBams {
    input:
      bams = ExtractBamsFromDir.bam_files
  }

  # Conditional QC execution
  if (run_QC || run_postqc) {
    call RunCallableLoci {
      input:
        indexed_bams = GetIndexedBams.indexed_bams
    }
    if (run_QC) {
      call QC {
        input:
          indexed_bams = GetIndexedBams.indexed_bams
      }
    }
  }

  # Step 4: Handle short variants, CNVs, and SVs
  if (run_snvs || run_indels || run_svs || run_cnvs || run_postqc) {
    call GetGzippedVcfs {
      input:
        vcfs = ExtractInputVcfFromDir.vcf_files
    }

    call GetGermlineVcfs {
      input:
        gzipped_vcfs = GetGzippedVcfs.gzipped_vcfs
    }

    if (run_postqc) {
      if (defined(optional_postqc_dirs["ptato_vcfs_dir"]) && defined(optional_postqc_dirs["walker_vcfs_dir"])) {
        call ExtractCombinedPtatoVcfFromDir {
          input:
            ptato_vcfs_dir = optional_postqc_dirs["ptato_vcfs_dir"]
        }

        call ExtractWalkerVcfFromDir {
          input:
            walker_vcfs_dir = optional_postqc_dirs["walker_vcfs_dir"]
        }

        call GetGzippedVcfs as GetGzippedWalkerVcfs {
          input:
            vcfs = ExtractWalkerVcfFromDir.walker_vcfs
        }

        call ExtractPtatoTableFromDir {
          input:
            ptato_vcfs_dir = optional_postqc_dirs["ptato_vcfs_dir"]
        }

        call PostPtatoQC {
          input:
            combined_vcfs = ExtractCombinedPtatoVcfFromDir.combined_vcfs,
            ptato_tables = ExtractPtatoTableFromDir.ptato_tables,
            walker_vcfs = GetGzippedWalkerVcfs.gzipped_vcfs,
            callable_loci_output = RunCallableLoci.callable_loci_output
        }
      } else {
        call ShortVariants {
          input:
            input_vcfs = GetGzippedVcfs.gzipped_vcfs,
            indexed_bams = GetIndexedBams.indexed_bams,
            germline_vcfs = GetGermlineVcfs.germline_vcfs
        }

        call PostPtatoQC {
          input:
            combined_vcfs = ShortVariants.output_vcfs,
            callable_loci_output = RunCallableLoci.callable_loci_output
        }
      }
    } else {
      if (run_snvs || run_indels) {
        call ShortVariants {
          input:
            input_vcfs = GetGzippedVcfs.gzipped_vcfs,
            indexed_bams = GetIndexedBams.indexed_bams,
            germline_vcfs = GetGermlineVcfs.germline_vcfs
        }
      }
    }

    if (run_svs || run_cnvs) {
      call CNVs {
        input:
          indexed_bams = GetIndexedBams.indexed_bams,
          germline_vcfs = GetGermlineVcfs.germline_vcfs
      }

      if (run_svs) {
        call SVs {
          input:
            indexed_bams = GetIndexedBams.indexed_bams,
            germline_vcfs = GetGermlineVcfs.germline_vcfs,
            cnv_files = CNVs.cnv_files
        }
      }
    }
  }

  output {
    Array[File] indexed_bams = GetIndexedBams.indexed_bams
    Array[File] gzipped_vcfs = GetGzippedVcfs.gzipped_vcfs
    Array[File] germline_vcfs = GetGermlineVcfs.germline_vcfs
  }
}
