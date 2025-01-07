version 1.0

import "../../NextflowModules/walker/2.2.0/walker.wdl" as Walker
import "../../NextflowModules/htslib/1.15/bgzip.wdl" as HtslibBgzip
import "../../NextflowModules/htslib/1.15/tabix.wdl" as HtslibTabix

workflow get_walker_vcfs {
  input {
    Array[Array[File]] input_somatic_vcfs
    Array[Array[File]] input_germline_vcfs
    Array[Array[File]] input_bams
    String out_dir
  }

  # Join germline VCFs with BAMs, then combine with somatic VCFs
  Array[Array[Array[File]]] input_walker = zip(input_germline_vcfs, input_bams).combine(input_somatic_vcfs)

  # Call Walker tool
  call Walker.walker {
    input:
      input_files = input_walker
  }

  # Call Bgzip on walker output VCFs
  call HtslibBgzip.bgzip {
    input:
      input_files = Walker.walker.out.map {
        Array[File] [
          it[0],  # donor_id
          it[1],  # sample_id
          it[2],  # walker_vcf
          sep('/', out_dir, "intermediate", "short_variants", "walker", ~{it[0]}, basename(it[3])),  # walker_bed
          sep('/', out_dir, "intermediate", "short_variants", "walker", ~{it[0]}, basename(it[4]))   # walker_txt
        ]
      }
  }

  # Call Tabix on bgzipped VCFs
  call HtslibTabix.tabix {
    input:
      input_files = HtslibBgzip.bgzip.out
  }

  # Copy and emit the final Walker VCFs and their index files
  Array[Array[File]] walker_vcfs = HtslibTabix.tabix.out.map {
    Array[File] [
      it[0],  # donor_id
      it[1],  # sample_id
      sep('/', out_dir, "intermediate", "short_variants", "walker", ~{it[0]}, basename(it[2])),  # vcf_gz
      sep('/', out_dir, "intermediate", "short_variants", "walker", ~{it[0]}, basename(it[3]))   # vcf_tbi
    ]
  }

  output {
    Array[Array[File]] walker_vcfs
  }
}
