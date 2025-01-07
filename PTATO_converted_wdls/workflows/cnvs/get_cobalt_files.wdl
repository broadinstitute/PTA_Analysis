version 1.0

workflow get_cobalt_files {
  input {
    Array[File] normal_bams
    Array[File] tumor_bams
    String out_dir
  }

  # Step 1: Combine normal and tumor BAM files
  scatter (i in range(length(normal_bams))) {
    call CountBamLinesApplication {
      input:
        normal_bam = normal_bams[i]
        tumor_bam = tumor_bams[i]
    }
  }

  # Step 2: Process and store cobalt ratio TSVs
  Array[Array[String]] cobalt_ratio_tsvs = CountBamLinesApplication.out.map {
    (donor_id, normal_sample_id, tumor_sample_id, cobalt_file, cobalt_ratio_tsv) =>
      String file_name = basename(cobalt_file)
      String tsv_name = basename(cobalt_ratio_tsv)
      File cobalt_file_copy = localizeFile(cobalt_file, "${out_dir}/intermediate/cnvs/cobalt/${donor_id}/${normal_sample_id}/${tumor_sample_id}/${file_name}")
      File cobalt_ratio_tsv_copy = localizeFile(cobalt_ratio_tsv, "${out_dir}/intermediate/cnvs/cobalt/${donor_id}/${normal_sample_id}/${tsv_name}")
      return [donor_id, normal_sample_id, tumor_sample_id, cobalt_ratio_tsv_copy]
  }

  output {
    Array[Array[String]] cobalt_ratio_tsvs
  }
}
