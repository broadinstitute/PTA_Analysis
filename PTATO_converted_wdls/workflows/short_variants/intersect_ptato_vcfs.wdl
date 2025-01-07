version 1.0

import "../../NextflowModules/bedtools/2.30.0/intersectPTATO.wdl" as BedtoolsIntersectPTATO
import "../../NextflowModules/htslib/1.15/bgzip.wdl" as HtslibBgzip
import "../../NextflowModules/htslib/1.15/tabix.wdl" as HtslibTabix

workflow intersect_ptato_vcfs {
  input {
    Array[Array[File]] input_vcfs
    Array[Array[File]] snvs_ptato_vcfs
    Array[Array[File]] indels_ptato_vcfs
    String out_dir
  }

  # Combine input VCFs for intersection
  scatter (intersect_vcf in input_vcfs) {
    Array[Array[File]] input_intersect_vcfs = zip(
      intersect_vcf,
      snvs_ptato_vcfs.zip(indels_ptato_vcfs).groupTuple(by: 0)
    )

    # Call the intersectPTATO tool
    call BedtoolsIntersectPTATO.intersectPTATO {
      input:
        input_files = input_intersect_vcfs
    }

    # Compress the intersected VCF using bgzip
    call HtslibBgzip.bgzip {
      input:
        input_files = BedtoolsIntersectPTATO.intersectPTATO.out
    }

    # Index the compressed VCF using tabix
    call HtslibTabix.tabix {
      input:
        input_files = HtslibBgzip.bgzip.out
    }

    # Prepare the final output
    Array[Array[File]] ptato_intersect_vcfs = HtslibTabix.tabix.out.map {
      Array[File] [
        it[0],  # donor_id
        it[1],  # sample_id
        sep('/', out_dir, "ptato_intersect_vcfs", ~{it[0]}, basename(it[2])),  # vcf_gz
        sep('/', out_dir, "ptato_intersect_vcfs", ~{it[0]}, basename(it[3]))   # vcf_tbi
      ]
    }
  }

  output {
    Array[Array[File]] ptato_intersect_vcfs
  }
}
