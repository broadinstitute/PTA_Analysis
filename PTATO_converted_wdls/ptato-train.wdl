version 1.0

workflow ptato-train {
  input {
    Array[String] bulk_names           # List of donor IDs (bulk names)
    String input_vcfs_dir              # Directory containing input VCFs
    Boolean run_snvs                   # Whether to run SNV analysis
    Boolean run_indels                 # Whether to run Indel analysis
    Boolean run_svs                    # Whether to run SV analysis
  }

  # Step 1: Generate donor IDs
  scatter (bulk_name in bulk_names) {
    call GenerateDonorIds {
      input:
        bulk_name = bulk_name
    }
  }

  # Step 2: Extract gzipped VCFs
  call ExtractInputVcfGzFromDir {
    input:
      input_vcfs_dir = input_vcfs_dir
  }

  # Step 3: Combine donor IDs with extracted VCFs
  call CombineInputVcfs {
    input:
      donor_ids = GenerateDonorIds.donor_ids,
      input_vcfs = ExtractInputVcfGzFromDir.vcf_files
  }

  # Step 4: Get gzipped VCFs
  call GetGzippedVcfs {
    input:
      input_vcfs = CombineInputVcfs.combined_vcfs
  }

  # Conditional Execution: Run Germline and Short Variant analysis
  if (run_snvs || run_indels || run_svs) {
    call GetGermlineVcfs {
      input:
        input_vcfs = GetGzippedVcfs.gzipped_vcfs
    }

    if (run_snvs || run_indels) {
      call ShortVariantsTrain {
        input:
          input_vcfs = GetGzippedVcfs.gzipped_vcfs,
          germline_vcfs = GetGermlineVcfs.germline_vcfs
      }
    }
  }
  
  output {
    Array[File] gzipped_vcfs = GetGzippedVcfs.gzipped_vcfs
    Array[File]? germline_vcfs = if (run_snvs || run_indels || run_svs) then GetGermlineVcfs.germline_vcfs else none
    Array[File]? trained_variants = if (run_snvs || run_indels) then ShortVariantsTrain.trained_variants else none
  }
}
task GenerateDonorIds {
  input {
    String bulk_name
  }
  command {
    echo "~{bulk_name}" > donor_id.txt
  }
  output {
    File donor_id = "donor_id.txt"
  }
}
task ExtractInputVcfGzFromDir {
  input {
    String input_vcfs_dir
  }
  command <<<
    find ~{input_vcfs_dir} -name "*.vcf.gz" > vcfs_list.txt
  >>>
  output {
    Array[File] vcf_files = read_lines("vcfs_list.txt")
  }
}
task CombineInputVcfs {
  input {
    Array[File] donor_ids
    Array[File] input_vcfs
  }
  command {
    paste <(echo ~{sep='\n' donor_ids}) <(echo ~{sep='\n' input_vcfs}) > combined_vcfs.txt
  }
  output {
    Array[File] combined_vcfs = read_lines("combined_vcfs.txt")
  }
}
task GetGzippedVcfs {
  input {
    Array[File] input_vcfs
  }
  command <<<
    gzip -c ~{sep=' ' input_vcfs} > gzipped_vcfs.txt
  >>>
  output {
    Array[File] gzipped_vcfs = read_lines("gzipped_vcfs.txt")
  }
}
task GetGermlineVcfs {
  input {
    Array[File] input_vcfs
  }
  command {
    # Example command for extracting germline VCFs
    echo "Extracting germline VCFs"
  }
  output {
    Array[File] germline_vcfs = ["germline1.vcf", "germline2.vcf"]
  }
}
task ShortVariantsTrain {
  input {
    Array[File] input_vcfs
    Array[File] germline_vcfs
  }
  command {
    echo "Training short variants with input VCFs and germline VCFs"
  }
  output {
    Array[File] trained_variants = ["variant1.vcf", "variant2.vcf"]
  }
}
