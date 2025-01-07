version 1.0

# Import external WDL modules required by the main workflow

import "./get_gzipped_vcfs.wdl" as GzipUtils
import "../NextflowModules/Utils/getFilesFromDir.wdl" as FileUtils
import "./short_variants/get_ab_tables.wdl" as AbTables
import "./short_variants/get_walker_vcfs.wdl" as WalkerUtils
import "./short_variants/get_somatic_vcfs.wdl" as SomaticUtils
import "./short_variants/get_context_beds.wdl" as ContextUtils
import "./short_variants/get_features_beds.wdl" as FeatureUtils
import "../NextflowModules/GATK/4.2.6.1/SplitVcfs.wdl" as GATKUtils
import "../NextflowModules/htslib/1.15/bgzip.wdl" as HtslibUtils
import "../NextflowModules/htslib/1.15/tabix.wdl" as TabixUtils

import "./short_variants/intersect_ptato_vcfs.wdl" as PtatoIntersect
import "./short_variants/merge_ptato_vcfs.wdl" as PtatoMerge

import "./snvs.wdl" as Snvs
import "./indels.wdl" as Indels

task extractSomaticVcf {
  input {
    String somatic_vcfs_dir
  }
  command {
    # Command to extract somatic VCFs
    ./extractSomaticVcfFromDir.sh ~{somatic_vcfs_dir}
  }
  output {
    Array[File] somatic_vcfs
  }
}

task getGzippedVcfs {
  input {
    Array[File] raw_somatic_vcfs
  }
  command {
    # Command to bgzip and tabix VCFs
    ./bgzip_tabix.sh ~{sep="," raw_somatic_vcfs}
  }
  output {
    Array[File] gzipped_vcfs
  }
}

task getSomaticVcfs {
  input {
    Array[File] input_vcfs
    Array[File] bams
  }
  command {
    # Command to generate somatic VCFs
    ./getSomaticVcfs.sh ~{sep="," input_vcfs} ~{sep="," bams}
  }
  output {
    Array[File] somatic_vcfs
  }
}

task combineSomaticVcfs {
  input {
    Array[File] somatic_vcfs
    Array[File] bams
  }
  command {
    # Command to combine somatic VCFs with BAMs
    ./combineSomaticVcfs.sh ~{sep="," somatic_vcfs} ~{sep="," bams}
  }
  output {
    Array[File] combined_somatic_vcfs
  }
}

workflow short_variants {
  input {
    Array[File] input_vcfs
    Array[File] bams
    Array[File] germline_vcfs
    String? context_beds_dir
    String? features_beds_dir
    Array[String]? closest_features
    Array[String]? intersect_features
  }

  # Get context beds
  if (defined(context_beds_dir)) {
    call ContextUtils.extractContextBedFromDir {
      input:
        context_beds_dir = select_first([context_beds_dir])
    }
  }
  if (!defined(context_beds_dir)) {
    call ContextUtils.get_context_beds {
      input:
        somatic_vcfs = somatic_vcfs
    }
  }
  Array[File] context_beds = select_first([extractContextBedFromDir.beds, get_context_beds.beds])

  # Get feature beds
  if (defined(features_beds_dir)) {
    call FeatureUtils.extractFeaturesBedFromDir {
      input:
        features_beds_dir = select_first([features_beds_dir])
    }
  }
  if (!defined(features_beds_dir)) {
    scatter (context_bed in context_beds) {
      call FeatureUtils.closest_feature {
        input:
          context_bed = context_bed,
          closest_features = select_first([closest_features])
      }
      
      call FeatureUtils.intersect_feature {
        input:
          context_bed = context_bed,
          intersect_features = select_first([intersect_features])
      }
    }
    
    call FeatureUtils.merge_features {
      input:
        context_beds = context_beds,
        closest_features = closest_feature.out,
        intersect_features = intersect_feature.out
    }
  }
  Array[File] features_beds = select_first([extractFeaturesBedFromDir.beds, merge_features.beds])

  # Split VCFs
  call GATKUtils.SplitVcfs {
    input:
      vcfs = somatic_vcfs
  }

  # Process SNVs and Indels
  if (params.run.snvs) {
    call Snvs.snvs {
      input:
        ab_tables = ab_tables,
        features_beds = features_beds,
        somatic_vcfs = SplitVcfs.snv_vcfs,
        walker_vcfs = walker_vcfs
    }
  }

  if (params.run.indels) {
    call Indels.indels {
      input:
        somatic_vcfs = SplitVcfs.indel_vcfs,
        context_beds = context_beds
    }
  }

  # Intersect and merge results
  call PtatoIntersect.intersect_ptato_vcfs {
    input:
      input_vcfs = input_vcfs,
      snvs_vcfs = select_first([snvs.ptato_vcfs, []]),
      indels_vcfs = select_first([indels.ptato_vcfs, []])
  }

  call PtatoMerge.merge_ptato_vcfs {
    input:
      intersect_vcfs = intersect_ptato_vcfs.out,
      snvs_vcfs = select_first([snvs.ptato_vcfs, []]),
      indels_vcfs = select_first([indels.ptato_vcfs, []])
  }

  output {
    Array[File] postqc_combined_input = merge_ptato_vcfs.combined_vcfs
  }
}

workflow short_variants_train {
  input {
    Array[File] input_vcfs
    Array[File] germline_vcfs
    String pta_vcfs_dir
    String nopta_vcfs_dir
    String? ab_tables_dir
    String? context_beds_dir
    String? features_beds_dir
    Array[String]? closest_features
    Array[String]? intersect_features
    Boolean run_snvs
    Boolean run_indels
  }

  # Extract VCFs
  call extractPtaVcfFromDir { input: vcfs_dir = pta_vcfs_dir }
  call extractNoptaVcfFromDir { input: vcfs_dir = nopta_vcfs_dir }
  
  Array[File] raw_input_vcfs = flatten([extractPtaVcfFromDir.vcfs, extractNoptaVcfFromDir.vcfs])
  
  call getGzippedVcfs { input: raw_vcfs = raw_input_vcfs }
  Array[File] somatic_vcfs = getGzippedVcfs.gzipped_vcfs

  # Get AB tables
  if (defined(ab_tables_dir)) {
    call extractAbTableFromDir { input: ab_tables_dir = select_first([ab_tables_dir]) }
  }
  if (!defined(ab_tables_dir)) {
    call get_ab_tables { 
      input:
        germline_vcfs = germline_vcfs,
        somatic_vcfs = somatic_vcfs
    }
  }
  Array[File] ab_tables = select_first([extractAbTableFromDir.tables, get_ab_tables.out])

    if ( params.optional.short_variants.features_beds_dir ) {
      features_beds = extractFeaturesBedFromDir( params.optional.short_variants.features_beds_dir )
    } else {
      closest_features = Channel.from( params.features.closest )
      intersect_features = Channel.from( params.features.intersect )

  # Process SNVs and Indels
  if (run_snvs) {
    call snvs_train {
      input:
        ab_tables = ab_tables,
        features_beds = features_beds,
        label_info = label_info
    }
  }

  if (run_indels) {
    call indels_train {
      input:
        ab_tables = ab_tables,
        features_beds = features_beds,
        label_info = label_info
    }
  }
}
