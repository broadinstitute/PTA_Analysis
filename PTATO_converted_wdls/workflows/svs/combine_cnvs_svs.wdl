version 1.0

import "../../NextflowModules/Utils/getFilesFromDir.wdl" as GetFilesUtils
import "../get_gzipped_vcfs.wdl" as GetGzippedVcfs
import "../../NextflowModules/Utils/filterGripss.wdl" as FilterGripss
import "../../NextflowModules/Utils/integrateSvFiles.wdl" as IntegrateSvFiles
import "../../NextflowModules/Utils/createSvPlots.wdl" as CreateSvPlots
import "../../NextflowModules/Utils/createCircosConfig.wdl" as CreateCircosConfig
import "../../NextflowModules/gridss-purple-linx/1.3.2/Circos.wdl" as Circos

workflow combine_cnvs_svs {

  input {
    Array[File] filtered_cnv_files
    Array[File] gripss_somatic_filtered_vcfs
    String? gripss_filtered_files_dir
    String? integrated_sv_files_dir
    String out_dir
  }

  call filterGripssFiles {
    input:
      filtered_cnv_files = filtered_cnv_files,
      gripss_somatic_filtered_vcfs = gripss_somatic_filtered_vcfs,
      gripss_filtered_files_dir = gripss_filtered_files_dir,
      out_dir = out_dir
  }

  call integrateSvFiles {
    input:
      filtered_cnv_files = filterGripssFiles.filtered_cnv_output,
      gripss_somatic_filtered_vcfs = filterGripssFiles.integrated_sv_output,
      integrated_sv_files_dir = integrated_sv_files_dir,
      out_dir = out_dir
  }

  output {
    Array[File] integrated_cnv_files = integrateSvFiles.integrated_cnv_files
    Array[File] integrated_sv_filtered_files = integrateSvFiles.integrated_sv_filtered_files
  }
}

task filterGripssFiles {
  input {
    Array[File] filtered_cnv_files
    Array[File] gripss_somatic_filtered_vcfs
    String? gripss_filtered_files_dir
    String out_dir
  }
  command {
    # Filtering logic for Gripss files goes here
    # Copy files to output directory
  }
  output {
    Array[File] filtered_cnv_output
  }
}

task integrateSvFiles {
  input {
    Array[File] filtered_cnv_files
    Array[File] gripss_somatic_filtered_vcfs
    String? integrated_sv_files_dir
    String out_dir
  }
  command {
    # Integration logic for SV files goes here
    # Copy files to output directory
  }
  output {
    Array[File] integrated_cnv_files
    Array[File] integrated_sv_filtered_files
  }
}

workflow combine_cnvs_svs {

  input {
    Array[File] filtered_cnv_files
    Array[File] gripss_somatic_filtered_vcfs
    String? gripss_filtered_files_dir
    String? integrated_sv_files_dir
    String out_dir
  }

  call filterGripssFiles {
    input:
      filtered_cnv_files = filtered_cnv_files,
      gripss_somatic_filtered_vcfs = gripss_somatic_filtered_vcfs,
      gripss_filtered_files_dir = gripss_filtered_files_dir,
      out_dir = out_dir
  }

  call integrateSvFiles {
    input:
      filtered_cnv_files = filterGripssFiles.filtered_cnv_output,
      gripss_somatic_filtered_vcfs = filterGripssFiles.integrated_sv_output,
      integrated_sv_files_dir = integrated_sv_files_dir,
      out_dir = out_dir
  }

  output {
    Array[File] integrated_cnv_files = integrateSvFiles.integrated_cnv_files
    Array[File] integrated_sv_filtered_files = integrateSvFiles.integrated_sv_filtered_files
  }
}

task filterGripssFiles {
  input {
    Array[File] filtered_cnv_files
    Array[File] gripss_somatic_filtered_vcfs
    String? gripss_filtered_files_dir
    String out_dir
  }
  command {
    # Filtering logic for Gripss files goes here
    # Copy files to output directory
  }
  output {
    Array[File] filtered_cnv_output
  }
}

task integrateSvFiles {
  input {
    Array[File] filtered_cnv_files
    Array[File] gripss_somatic_filtered_vcfs
    String? integrated_sv_files_dir
    String out_dir
  }
  command {
    # Integration logic for SV files goes here
    # Copy files to output directory
  }
  output {
    Array[File] integrated_cnv_files
    Array[File] integrated_sv_filtered_files
  }
}
