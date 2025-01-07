version 1.0

import "../../NextflowModules/Utils/mergePtatoVcfs.wdl" as MergePtato
import "../../NextflowModules/htslib/1.15/bgzip.wdl" as HtslibBgzip
import "../../NextflowModules/htslib/1.15/tabix.wdl" as HtslibTabix

workflow merge_ptato_vcfs {
  input {
    Array[Array[File]] ptato_intersect_vcfs
    Array[Array[File]] snvs_ptato_vcfs
    Array[Array[File]] indels_ptato_vcfs
    String out_dir
  }

  # Combine input VCFs for merging
  scatter (intersect_vcf in ptato_intersect_vcfs) {
    Array[Array[File]] input_merge_vcfs = zip(
      intersect_vcf,
      snvs_ptato_vcfs.zip(indels_ptato_vcfs).groupTuple(by: 0)
    )
    
    # Call the mergePtatoVcfs tool
    call MergePtato.mergePtatoVcfs {
      input:
        input_files = input_merge_vcfs
    }
    
    # Compress the merged VCF using bgzip
    call HtslibBgzip.bgzip {
      input:
        input_files = MergePtato.mergePtatoVcfs.out
    }
    
    # Index the compressed VCF using tabix
    call HtslibTabix.tabix {
      input:
        input_files = HtslibBgzip.bgzip.out
    }
    
    # Copy and organize the output files
    Array[Array[File]] ptato_merged_vcfs = HtslibTabix.tabix.out.map {
      Array[File] [
        it[0],  # donor_id
        it[1],  # sample_id
        sep('/', out_dir, "ptato_vcfs", ~{it[0]}, basename(it[2])),  # vcf_gz
        sep('/', out_dir, "ptato_vcfs", ~{it[0]}, basename(it[3]))   # vcf_tbi
      ]
    }
  }

  output {
    Array[Array[File]] ptato_merged_vcfs
  }
}
