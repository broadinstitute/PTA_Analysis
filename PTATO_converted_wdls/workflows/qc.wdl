version 1.0

import "../NextflowModules/Utils/getFilesFromDir.wdl" as Utils
import "../NextflowModules/GATK/4.2.6.1/CollectWGSMetrics.wdl" as GATK_WGS
import "../NextflowModules/GATK/4.2.6.1/CollectAlignmentSummaryMetrics.wdl" as GATK_ASM
import "../NextflowModules/GATK/3.8.1/CallableLoci.wdl" as GATK_CallableLoci
import "../NextflowModules/Utils/QCreport.wdl" as QC
import "../NextflowModules/Utils/postQCreport.wdl" as PostQC

workflow qc {
  input {
    Array[File] bams
    String out_dir
    String? alignment_summary_metrics_dir
    String? wgs_metrics_dir
  }

  # Step 1: Collect alignment summary metrics or extract if provided
  if (defined(alignment_summary_metrics_dir)) {
    call Utils.extractAlignmentSummaryMetricsFromDir {
      input:
        dir = alignment_summary_metrics_dir
    }
    Array[Array[File]] alignment_summary_metrics_files = Utils.extractAlignmentSummaryMetricsFromDir.out
  } else {
    call GATK_ASM.CollectAlignmentSummaryMetrics {
      input:
        bams = bams
    }
    Array[Array[File]] alignment_summary_metrics_files = GATK_ASM.CollectAlignmentSummaryMetrics.out.map { 
      Array[File] [
        sep('/', out_dir, "QC", "alignment_summary_metrics", ~{it[0]}, basename(it[2]))
      ]
    }
  }

  # Step 2: Collect WGS metrics or extract if provided
  if (defined(wgs_metrics_dir)) {
    call Utils.extractWGSMetricsFromDir {
      input:
        dir = wgs_metrics_dir
    }
    Array[Array[File]] wgs_metrics_files = Utils.extractWGSMetricsFromDir.out
  } else {
    call GATK_WGS.CollectWGSMetrics {
      input:
        bams = bams
    }
    Array[Array[File]] wgs_metrics_files = GATK_WGS.CollectWGSMetrics.out.map { 
      Array[File] [
        sep('/', out_dir, "QC", "wgs_metrics", ~{it[0]}, basename(it[2]))
      ]
    }
  }

  # Step 3: Combine input for QC report
  Array[Array[File]] input_qc_report = zip(alignment_summary_metrics_files, wgs_metrics_files)

  call QC.QCreport {
    input:
      input_qc = input_qc_report
  }

  Array[Array[File]] qc_reports = QC.QCreport.out.map { 
    Array[File] [
      sep('/', out_dir, "QC", "reports", ~{it[0]}, basename(it[1])),
      sep('/', out_dir, "QC", "reports", ~{it[0]}, basename(it[2]))
    ]
  }

  output {
    Array[Array[File]] qc_reports
  }
}

workflow post_ptato_qc {
  input {
    Array[Array[File]] postqc_combined_input
    Array[Array[File]] callableloci_files
    String out_dir
  }

  # Step 1: Rename sample IDs for PostQC
  Array[Array[File]] renamed_postqc_input = postqc_combined_input.map { 
    Array[File] [
      it[0],
      it[1].replaceAll("_.*", ""),
      it[2],
      it[3],
      it[4],
      it[5],
      it[6],
      it[7],
      it[8]
    ]
  }

  Array[Array[File]] renamed_callableloci_files = callableloci_files.map { 
    Array[File] [
      it[0],
      it[1],
      it[4]
    ]
  }

  # Step 2: Combine input for PostQC report
  Array[Array[File]] input_PostPTATO_qc = zip(renamed_postqc_input, renamed_callableloci_files)

  call PostQC.postQCreport {
    input:
      input_postqc = input_PostPTATO_qc
  }

  Array[Array[File]] postQC_report = PostQC.postQCreport.out.map { 
    Array[File] [
      sep('/', out_dir, "QC", "reports", ~{it[0]}, basename(it[2])),
      sep('/', out_dir, "QC", "reports", ~{it[0]}, basename(it[3]))
    ]
  }

  output {
    Array[Array[File]] postQC_report
  }
}
