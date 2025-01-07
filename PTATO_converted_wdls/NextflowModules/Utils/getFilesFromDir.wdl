version 1.0
workflow GenomicFileExtractionWorkflow {
  input {
    String nopta_dir
    String pta_dir
    String input_dir
    String germline_vcfs_dir
    String somatic_dir
    String bams_dir
    String walker_dir
    String ab_dir
    String context_bed_dir
    String features_bed_dir
    String phased_dir
    String excludeindellist_filtered_vcfs_dir
    String snv_rf_tables_dir
    String indel_rf_tables_dir
    String ptato_vcfs_dir
    String ptato_table_dir
    String wgs_metrics_dir
    String alignment_summary_metrics_dir
    String autosomal_callable_loci_dir
    String callable_loci_bed_dir
    String gridss_driver_vcfs_dir
    String gridss_unfiltered_vcfs_dir
    String gripss_somatic_filtered_vcfs_dir
    String cobalt_ratio_tsv_dir
    String cobalt_filtered_readcounts_dir
    String baf_filtered_files_dir
    String baf_binned_files_dir
    String baf_segments_files_dir
    String gripss_filtered_files_dir
    String integrated_sv_files_dir
  }

  call ExtractNoptaVcf {
    input:
      directory = nopta_dir
  }

  call ExtractPtaVcf {
    input:
      directory = pta_dir
  }

  call ExtractInputVcf {
    input:
      directory = input_dir
  }

  call ExtractGermlineVcf {
    input:
      directory = germline_vcfs_dir
  }

  call ExtractSomaticVcf {
    input:
      directory = somatic_dir
  }

  call ExtractBams {
    input:
      directory = bams_dir
  }

  call ExtractWalkerVcf {
    input:
      directory = walker_dir
  }

  call ExtractAbTable {
    input:
      directory = ab_dir
  }

  call ExtractContextBed {
    input:
      directory = context_bed_dir
  }

  call ExtractFeaturesBed {
    input:
      directory = features_bed_dir
  }

  call ExtractPhasedVcfGz {
    input:
      directory = phased_dir
  }

  call ExtractExcludeIndelListFilteredVcf {
    input:
      directory = excludeindellist_filtered_vcfs_dir
  }

  call ExtractSnvRfTableRds {
    input:
      directory = snv_rf_tables_dir
  }

  call ExtractIndelRfTableRds {
    input:
      directory = indel_rf_tables_dir
  }

  call ExtractPtatoVcf {
    input:
      directory = ptato_vcfs_dir
  }

  call ExtractPtatoTable {
    input:
      directory = ptato_table_dir
  }

  call ExtractCombinedPtatoVcf {
    input:
      directory = ptato_vcfs_dir
  }

  call ExtractWGSMetrics {
    input:
      directory = wgs_metrics_dir
  }

  call ExtractAlignmentSummaryMetrics {
    input:
      directory = alignment_summary_metrics_dir
  }

  call ExtractAutosomalCallableLoci {
    input:
      directory = autosomal_callable_loci_dir
  }

  call ExtractCallableLociBed {
    input:
      directory = callable_loci_bed_dir
  }

  call ExtractGridssDriverVcf {
    input:
      directory = gridss_driver_vcfs_dir
  }

  call ExtractGridssUnfilteredVcf {
    input:
      directory = gridss_unfiltered_vcfs_dir
  }

  call ExtractGripssSomaticFilteredVcf {
    input:
      directory = gripss_somatic_filtered_vcfs_dir
  }

  call ExtractCobaltRatioTsv {
    input:
      directory = cobalt_ratio_tsv_dir
  }

  call ExtractCobaltFilteredReadCounts {
    input:
      directory = cobalt_filtered_readcounts_dir
  }

  call ExtractCobaltFilteredReadCountsSegments {
    input:
      directory = cobalt_filtered_readcounts_dir
  }

  call ExtractBafFilteredFiles {
    input:
      directory = baf_filtered_files_dir
  }

  call ExtractBafBinnedFiles {
    input:
      directory = baf_binned_files_dir
  }

  call ExtractBafSegmentsFiles {
    input:
      directory = baf_segments_files_dir
  }

  call ExtractGripssFilteredFiles {
    input:
      directory = gripss_filtered_files_dir
  }

  call ExtractIntegratedSvFiles {
    input:
      directory = integrated_sv_files_dir
  }

  output {
    Array[Array[String]] nopta_vcf_files = ExtractNoptaVcf.extracted_files
    Array[Array[String]] pta_vcf_files = ExtractPtaVcf.extracted_files
    Array[Array[String]] input_vcf_files = ExtractInputVcf.extracted_files
    Array[Array[String]] germline_vcf_files = ExtractGermlineVcf.extracted_files
    Array[Array[String]] somatic_vcf_files = ExtractSomaticVcf.extracted_files
    Array[Array[String]] bam_files = ExtractBams.extracted_files
    Array[Array[String]] walker_vcf_files = ExtractWalkerVcf.extracted_files
    Array[Array[String]] ab_table_files = ExtractAbTable.extracted_files
    Array[Array[String]] context_bed_files = ExtractContextBed.extracted_files
    Array[Array[String]] features_bed_files = ExtractFeaturesBed.extracted_files
    Array[Array[String]] phased_vcf_files = ExtractPhasedVcfGz.extracted_files
    Array[Array[String]] exclude_indel_list_filtered_vcf_files = ExtractExcludeIndelListFilteredVcf.extracted_files
    Array[Array[String]] snv_rf_table_rds_files = ExtractSnvRfTableRds.extracted_files
    Array[Array[String]] indel_rf_table_rds_files = ExtractIndelRfTableRds.extracted_files
    Array[Array[String]] ptato_vcf_files = ExtractPtatoVcf.extracted_files
    Array[Array[String]] ptato_table_files = ExtractPtatoTable.extracted_files
    Array[Array[String]] combined_ptato_vcf_files = ExtractCombinedPtatoVcf.extracted_files
    Array[Array[String]] wgs_metrics_files = ExtractWGSMetrics.extracted_files
    Array[Array[String]] alignment_summary_metrics_files = ExtractAlignmentSummaryMetrics.extracted_files
    Array[Array[String]] autosomal_callable_loci_files = ExtractAutosomalCallableLoci.extracted_files
    Array[Array[String]] callable_loci_bed_files = ExtractCallableLociBed.extracted_files
    Array[Array[String]] gridss_driver_vcf_files = ExtractGridssDriverVcf.extracted_files
    Array[Array[String]] gridss_unfiltered_vcf_files = ExtractGridssUnfilteredVcf.extracted_files
    Array[Array[String]] gripss_somatic_filtered_vcf_files = ExtractGripssSomaticFilteredVcf.extracted_files
    Array[Array[String]] cobalt_ratio_tsv_files = ExtractCobaltRatioTsv.extracted_files
    Array[Array[String]] cobalt_filtered_readcounts_files = ExtractCobaltFilteredReadCounts.extracted_files
    Array[Array[String]] cobalt_filtered_readcounts_segments_files = ExtractCobaltFilteredReadCountsSegments.extracted_files
    Array[Array[String]] baf_filtered_files = ExtractBafFilteredFiles.extracted_files
    Array[Array[String]] baf_binned_files = ExtractBafBinnedFiles.extracted_files
    Array[Array[String]] baf_segments_files = ExtractBafSegmentsFiles.extracted_files
    Array[Array[String]] gripss_filtered_files = ExtractGripssFilteredFiles.extracted_files
    Array[Array[String]] integrated_sv_files = ExtractIntegratedSvFiles.extracted_files
  }
}

task ExtractNoptaVcf {
  input {
    String directory
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    find ~{directory} -type f \( -name "*.vcf" -o -name "*.vcf.gz" \) > nopta_vcf_files.txt

    # Check if corresponding .tbi files exist
    awk '{print $0".tbi"}' nopta_vcf_files.txt | xargs -I {} sh -c 'if [ ! -f {} ]; then echo "{} missing"; fi'
  >>>

  output {
    Array[Array[String]] extracted_files = read_tsv("nopta_vcf_files.txt")
  }

  runtime {
    docker: "ubuntu:latest"
    cpu: 2
    memory: "4 GB"
  }
}

task ExtractPtaVcf {
  input {
    String directory
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    find ~{directory} -type f \( -name "*.vcf" -o -name "*.vcf.gz" \) > pta_vcf_files.txt

    # Check if corresponding .tbi files exist
    awk '{print $0".tbi"}' pta_vcf_files.txt | xargs -I {} sh -c 'if [ ! -f {} ]; then echo "{} missing"; fi'
  >>>

  output {
    Array[Array[String]] extracted_files = read_tsv("pta_vcf_files.txt")
  }

  runtime {
    docker: "ubuntu:latest"
    cpu: 2
    memory: "4 GB"
  }
}

########################
workflow GenomicFileExtractionWorkflow {
  input {
    String nopta_dir
    String pta_dir
    String input_dir
    String germline_vcfs_dir
  }

  call ExtractNoptaVcf {
    input:
      directory = nopta_dir
  }

  call ExtractPtaVcf {
    input:
      directory = pta_dir
  }

  call ExtractInputVcf {
    input:
      directory = input_dir
  }

  call ExtractGermlineVcf {
    input:
      directory = germline_vcfs_dir
  }

  output {
    Array[Array[String]] nopta_vcf_files = ExtractNoptaVcf.extracted_files
    Array[Array[String]] pta_vcf_files = ExtractPtaVcf.extracted_files
    Array[Array[String]] input_vcf_files = ExtractInputVcf.extracted_files
    Array[Array[String]] germline_vcf_files = ExtractGermlineVcf.extracted_files
  }
}


task ExtractInputVcf {
  input {
    String directory
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    find ~{directory} -type f \( -name "*.vcf" -o -name "*.vcf.gz" \) > input_vcf_files.txt

    # Check if corresponding .tbi files exist
    awk '{print $0".tbi"}' input_vcf_files.txt | xargs -I {} sh -c 'if [ ! -f {} ]; then echo "{} missing"; fi'
  >>>

  output {
    Array[Array[String]] extracted_files = read_tsv("input_vcf_files.txt")
  }

  runtime {
    docker: "ubuntu:latest"
    cpu: 2
    memory: "4 GB"
  }
}




task ExtractGermlineVcf {
  input {
    String directory
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    find ~{directory} -type f \( -name "*.vcf" -o -name "*.vcf.gz" \) > germline_vcf_files.txt

    # Check if corresponding .tbi files exist
    awk '{print $0".tbi"}' germline_vcf_files.txt | xargs -I {} sh -c 'if [ ! -f {} ]; then echo "{} missing"; fi'
  >>>

  output {
    Array[Array[String]] extracted_files = read_tsv("germline_vcf_files.txt")
  }

  runtime {
    docker: "ubuntu:latest"
    cpu: 2
    memory: "4 GB"
  }
}

task ExtractSomaticVcf {
  input {
    String directory
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    find ~{directory} -type f \( -name "*.vcf" -o -name "*.vcf.gz" \) > somatic_vcf_files.txt

    # Check if corresponding .tbi files exist
    awk '{print $0".tbi"}' somatic_vcf_files.txt | xargs -I {} sh -c 'if [ ! -f {} ]; then echo "{} missing"; fi'
  >>>

  output {
    Array[Array[String]] extracted_files = read_tsv("somatic_vcf_files.txt")
  }

  runtime {
    docker: "ubuntu:latest"
    cpu: 2
    memory: "4 GB"
  }
}

task ExtractBams {
  input {
    String directory
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    find ~{directory} -type f -name "*.bam" > bam_files.txt

    # Check if corresponding .bai files exist
    awk '{print $0".bai"}' bam_files.txt | xargs -I {} sh -c 'if [ -f {} ]; then echo "{} exists"; else echo "{} missing"; fi' > bam_indices_check.txt
  >>>

  output {
    Array[Array[String]] extracted_files = read_tsv("bam_files.txt")
  }

  runtime {
    docker: "ubuntu:latest"
    cpu: 2
    memory: "4 GB"
  }
}

task ExtractWalkerVcf {
  input {
    String directory
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    find ~{directory} -type f \( -name "*.vcf" -o -name "*.vcf.gz" \) > walker_vcf_files.txt

    # Check if corresponding .tbi files exist
    awk '{print $0".tbi"}' walker_vcf_files.txt | xargs -I {} sh -c 'if [ ! -f {} ]; then echo "{} missing"; fi'
  >>>

  output {
    Array[Array[String]] extracted_files = read_tsv("walker_vcf_files.txt")
  }

  runtime {
    docker: "ubuntu:latest"
    cpu: 2
    memory: "4 GB"
  }
}

task ExtractAbTable {
  input {
    String directory
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    find ~{directory} -type f -name "*.abtable.txt" > ab_table_files.txt
  >>>

  output {
    Array[Array[String]] extracted_files = read_tsv("ab_table_files.txt")
  }

  runtime {
    docker: "ubuntu:latest"
    cpu: 2
    memory: "4 GB"
  }
}

task ExtractContextBed {
  input {
    String directory
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    find ~{directory} -type f -name "*.bed" > context_bed_files.txt
  >>>

  output {
    Array[Array[String]] extracted_files = read_tsv("context_bed_files.txt")
  }

  runtime {
    docker: "ubuntu:latest"
    cpu: 2
    memory: "4 GB"
  }
}

task ExtractFeaturesBed {
  input {
    String directory
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    find ~{directory} -type f -name "*.bed" > features_bed_files.txt
  >>>

  output {
    Array[Array[String]] extracted_files = read_tsv("features_bed_files.txt")
  }

  runtime {
    docker: "ubuntu:latest"
    cpu: 2
    memory: "4 GB"
  }
}

task ExtractPhasedVcfGz {
  input {
    String directory
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    find ~{directory} -type f -name "*.phased.vcf.gz" > phased_vcf_files.txt

    # Check if corresponding .tbi files exist
    awk '{print $0".tbi"}' phased_vcf_files.txt | xargs -I {} sh -c 'if [ ! -f {} ]; then echo "{} missing"; fi'
  >>>

  output {
    Array[Array[String]] extracted_files = read_tsv("phased_vcf_files.txt")
  }

  runtime {
    docker: "ubuntu:latest"
    cpu: 2
    memory: "4 GB"
  }
}

task ExtractExcludeIndelListFilteredVcf {
  input {
    String directory
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    find ~{directory} -type f \( -name "*.vcf" -o -name "*.vcf.gz" \) > exclude_indel_list_filtered_vcfs.txt

    # Check if corresponding .tbi files exist
    awk '{print $0".tbi"}' exclude_indel_list_filtered_vcfs.txt | xargs -I {} sh -c 'if [ ! -f {} ]; then echo "{} missing"; fi'
  >>>

  output {
    Array[Array[String]] extracted_files = read_tsv("exclude_indel_list_filtered_vcfs.txt")
  }

  runtime {
    docker: "ubuntu:latest"
    cpu: 2
    memory: "4 GB"
  }
}

task ExtractSnvRfTableRds {
  input {
    String directory
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    find ~{directory} -type f -name "*.rftable.rds" > snv_rf_tables.txt
  >>>

  output {
    Array[Array[String]] extracted_files = read_tsv("snv_rf_tables.txt")
  }

  runtime {
    docker: "ubuntu:latest"
    cpu: 2
    memory: "4 GB"
  }
}
task ExtractIndelRfTableRds {
  input {
    String directory
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    find ~{directory} -type f -name "*.rftable.rds" > indel_rf_tables.txt
  >>>

  output {
    Array[Array[String]] extracted_files = read_tsv("indel_rf_tables.txt")
  }

  runtime {
    docker: "ubuntu:latest"
    cpu: 2
    memory: "4 GB"
  }
}

task ExtractPtatoVcf {
  input {
    String directory
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    find ~{directory} -type f \( -name "*.vcf" -o -name "*.vcf.gz" \) > ptato_vcf_files.txt

    # Check if corresponding .tbi files exist
    awk '{print $0".tbi"}' ptato_vcf_files.txt | xargs -I {} sh -c 'if [ ! -f {} ]; then echo "{} missing"; fi'
  >>>

  output {
    Array[Array[String]] extracted_files = read_tsv("ptato_vcf_files.txt")
  }

  runtime {
    docker: "ubuntu:latest"
    cpu: 2
    memory: "4 GB"
  }
}
task ExtractPtatoTable {
  input {
    String directory
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    find ~{directory} -type f -name "*.txt" > ptato_table_files.txt
  >>>

  output {
    Array[Array[String]] extracted_files = read_tsv("ptato_table_files.txt")
  }

  runtime {
    docker: "ubuntu:latest"
    cpu: 2
    memory: "4 GB"
  }
}
task ExtractCombinedPtatoVcf {
  input {
    String directory
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    find ~{directory} -type f \( -name "*.vcf" -o -name "*.vcf.gz" \) > combined_ptato_vcf_files.txt

    # Check if corresponding .tbi files exist
    awk '{print $0".tbi"}' combined_ptato_vcf_files.txt | xargs -I {} sh -c 'if [ ! -f {} ]; then echo "{} missing"; fi'
  >>>

  output {
    Array[Array[String]] extracted_files = read_tsv("combined_ptato_vcf_files.txt")
  }

  runtime {
    docker: "ubuntu:latest"
    cpu: 2
    memory: "4 GB"
  }
}
task ExtractWGSMetrics {
  input {
    String directory
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    find ~{directory} -type f -name "*.wgs_metrics*" > wgs_metrics_files.txt
  >>>

  output {
    Array[Array[String]] extracted_files = read_tsv("wgs_metrics_files.txt")
  }

  runtime {
    docker: "ubuntu:latest"
    cpu: 2
    memory: "4 GB"
  }
}

task ExtractAlignmentSummaryMetrics {
  input {
    String directory
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    find ~{directory} -type f -name "*.alignment_summary_metrics*" > alignment_summary_metrics_files.txt
  >>>

  output {
    Array[Array[String]] extracted_files = read_tsv("alignment_summary_metrics_files.txt")
  }

  runtime {
    docker: "ubuntu:latest"
    cpu: 2
    memory: "4 GB"
  }
}

task ExtractAutosomalCallableLoci {
  input {
    String directory
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    find ~{directory} -type f -name "*.callableloci.autosomal.txt" > autosomal_callable_loci_files.txt
  >>>

  output {
    Array[Array[String]] extracted_files = read_tsv("autosomal_callable_loci_files.txt")
  }

  runtime {
    docker: "ubuntu:latest"
    cpu: 2
    memory: "4 GB"
  }
}

task ExtractCallableLociBed {
  input {
    String directory
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    find ~{directory} -type f -name "*.callableloci.bed" > callable_loci_bed_files.txt

    # Check for corresponding .txt files
    awk '{print $0".txt"}' callable_loci_bed_files.txt | xargs -I {} sh -c 'if [ ! -f {} ]; then echo "{} missing"; fi'
  >>>

  output {
    Array[Array[String]] extracted_files = read_tsv("callable_loci_bed_files.txt")
  }

  runtime {
    docker: "ubuntu:latest"
    cpu: 2
    memory: "4 GB"
  }
}

task ExtractGridssDriverVcf {
  input {
    String directory
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    find ~{directory} -type f \( -name "*.gridss.driver.vcf" -o -name "*.gridss.driver.vcf.gz" \) > gridss_driver_vcf_files.txt

    # Check if corresponding .tbi files exist
    awk '{print $0".tbi"}' gridss_driver_vcf_files.txt | xargs -I {} sh -c 'if [ ! -f {} ]; then echo "{} missing"; fi'
  >>>

  output {
    Array[Array[String]] extracted_files = read_tsv("gridss_driver_vcf_files.txt")
  }

  runtime {
    docker: "ubuntu:latest"
    cpu: 2
    memory: "4 GB"
  }
}

task ExtractGridssUnfilteredVcf {
  input {
    String directory
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    find ~{directory} -type f \( -name "*.gridss.unfiltered.vcf" -o -name "*.gridss.unfiltered.vcf.gz" \) > gridss_unfiltered_vcf_files.txt

    # Check if corresponding .tbi files exist
    awk '{print $0".tbi"}' gridss_unfiltered_vcf_files.txt | xargs -I {} sh -c 'if [ ! -f {} ]; then echo "{} missing"; fi'
  >>>

  output {
    Array[Array[String]] extracted_files = read_tsv("gridss_unfiltered_vcf_files.txt")
  }

  runtime {
    docker: "ubuntu:latest"
    cpu: 2
    memory: "4 GB"
  }
}

task ExtractGripssSomaticFilteredVcf {
  input {
    String directory
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    find ~{directory} -type f \( -name "*.gripss.somatic.filtered.vcf" -o -name "*.gripss.somatic.filtered.vcf.gz" \) > gripss_somatic_filtered_vcfs.txt

    # Check if corresponding .tbi files exist
    awk '{print $0".tbi"}' gripss_somatic_filtered_vcfs.txt | xargs -I {} sh -c 'if [ ! -f {} ]; then echo "{} missing"; fi'
  >>>

  output {
    Array[Array[String]] extracted_files = read_tsv("gripss_somatic_filtered_vcfs.txt")
  }

  runtime {
    docker: "ubuntu:latest"
    cpu: 2
    memory: "4 GB"
  }
}
task ExtractCobaltRatioTsv {
  input {
    String directory
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    find ~{directory} -type f -name "*.cobalt.ratio.tsv" > cobalt_ratio_files.txt
  >>>

  output {
    Array[Array[String]] extracted_files = read_tsv("cobalt_ratio_files.txt")
  }

  runtime {
    docker: "ubuntu:latest"
    cpu: 2
    memory: "4 GB"
  }
}

task ExtractCobaltFilteredReadCounts {
  input {
    String directory
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    find ~{directory} -type f -name "*.readcounts.*.txt" > cobalt_filtered_readcounts.txt
  >>>

  output {
    Array[Array[String]] extracted_files = read_tsv("cobalt_filtered_readcounts.txt")
  }

  runtime {
    docker: "ubuntu:latest"
    cpu: 2
    memory: "4 GB"
  }
}

task ExtractCobaltFilteredReadCountsSegments {
  input {
    String directory
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    find ~{directory} -type f -name "*.readcounts.segments.*" > cobalt_filtered_readcounts_segments.txt
  >>>

  output {
    Array[Array[String]] extracted_files = read_tsv("cobalt_filtered_readcounts_segments.txt")
  }

  runtime {
    docker: "ubuntu:latest"
    cpu: 2
    memory: "4 GB"
  }
}

task ExtractBafFilteredFiles {
  input {
    String directory
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    find ~{directory} -type f -name "*.baf.filtered.txt" > baf_filtered_files.txt
  >>>

  output {
    Array[Array[String]] extracted_files = read_tsv("baf_filtered_files.txt")
  }

  runtime {
    docker: "ubuntu:latest"
    cpu: 2
    memory: "4 GB"
  }
}
task ExtractBafBinnedFiles {
  input {
    String directory
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    find ~{directory} -type f -name "*.baf.binned*.txt" > baf_binned_files.txt
  >>>

  output {
    Array[Array[String]] extracted_files = read_tsv("baf_binned_files.txt")
  }

  runtime {
    docker: "ubuntu:latest"
    cpu: 2
    memory: "4 GB"
  }
}

task ExtractBafSegmentsFiles {
  input {
    String directory
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    find ~{directory} -type f -name "*.baf.segments.*" > baf_segments_files.txt
  >>>

  output {
    Array[Array[String]] extracted_files = read_tsv("baf_segments_files.txt")
  }

  runtime {
    docker: "ubuntu:latest"
    cpu: 2
    memory: "4 GB"
  }
}

task ExtractGripssFilteredFiles {
  input {
    String directory
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    find ~{directory} -type f \( -name "*.svs.postfilter.vcf" -o -name "*.svs.postfilter.vcf.gz" \) > gripss_filtered_files.txt
  >>>

  output {
    Array[Array[String]] extracted_files = read_tsv("gripss_filtered_files.txt")
  }

  runtime {
    docker: "ubuntu:latest"
    cpu: 2
    memory: "4 GB"
  }
}

task ExtractIntegratedSvFiles {
  input {
    String directory
  }

  command <<<
    set -euo pipefail
    HOSTNAME=$(hostname)
    echo "Running on host: ${HOSTNAME}"

    find ~{directory} -type f -name "*.integrated.*" > integrated_sv_files.txt
  >>>

  output {
    Array[Array[String]] extracted_files = read_tsv("integrated_sv_files.txt")
  }

  runtime {
    docker: "ubuntu:latest"
    cpu: 2
    memory: "4 GB"
  }
}

