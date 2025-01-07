version 1.0

import "../NextflowModules/htslib/1.15/bgzip.wdl" as bgzip_tasks
import "../NextflowModules/htslib/1.15/tabix.wdl" as tabix_tasks

workflow get_gzipped_vcfs {
    input {
        Array[File] input_vcfs
    }

    # Filter for uncompressed VCFs
    Array[File] input_unzipped_vcfs = select_all(input_vcfs[where vcf.basename.endsWith(".vcf")])
    
    # Filter for gzipped VCFs without index
    Array[File] input_unindexed_vcfs = select_all(input_vcfs[where vcf.basename.endsWith(".vcf.gz") && !has_index])
    
    # Filter for already gzipped and indexed VCFs
    Array[Pair[File, File]] gzipped_vcf_files = select_all(zip(input_vcfs, indices)[where vcf.basename.endsWith(".vcf.gz") && has_index])

    # Compress unzipped VCFs
    scatter (vcf in input_unzipped_vcfs) {
        call bgzip_tasks.bgzip { input: vcf = vcf }
    }

    # Combine newly compressed VCFs with previously unindexed ones
    Array[File] all_unindexed_vcfs = flatten([input_unindexed_vcfs, bgzip.compressed_vcf])

    # Create indices for all unindexed VCFs
    scatter (vcf in all_unindexed_vcfs) {
        call tabix_tasks.tabix { input: vcf_file = vcf }
    }

    # Combine all gzipped and indexed VCFs
    Array[Pair[File, File]] final_gzipped_vcfs = flatten([gzipped_vcf_files, zip(all_unindexed_vcfs, tabix.index)])

    output {
        Array[Pair[File, File]] gzipped_vcf_files = final_gzipped_vcfs
    }
}
