version 1.0

workflow GenomicsDBImportWorkflow {
  meta {
    author: "Shadi Zaheri"
    description: "Imports multiple GVCFs into a GenomicsDB workspace for joint genotyping."
  }

  parameter_meta {
    reference_fasta: "Reference genome FASTA."
    reference_fasta_index: "FASTA index (.fai)."
    reference_dict: "Sequence dictionary (.dict)."
    gvcf_files: "List of per-sample GVCF files."
    sample_names: "List of sample names matching the GVCF files."
    genomicsdb_workspace: "Name of the GenomicsDB workspace directory to create."
    interval_file: "Optional: Interval list for import."
    preemptible_tries: "Number of preemptible retries."
    memory_gb: "Memory per task (GiB)."
    disk_gb: "Disk per task (GiB)."
    cpu: "CPU cores per task."
    gatk_docker: "GATK Docker image (default: broadinstitute/gatk:4.6.1.0)."
  }

  input {
    File reference_fasta
    File reference_fasta_index
    File reference_dict
    Array[File] gvcf_files
    Array[String] sample_names
    String genomicsdb_workspace
    File? interval_file
    Int preemptible_tries = 1
    Int memory_gb = 32
    Int disk_gb = 200
    Int cpu = 4
    String gatk_docker = "broadinstitute/gatk:4.6.1.0"
  }

  call GenomicsDBImportTask {
    input:
      reference_fasta = reference_fasta,
      reference_fasta_index = reference_fasta_index,
      reference_dict = reference_dict,
      gvcf_files = gvcf_files,
      sample_names = sample_names,
      genomicsdb_workspace = genomicsdb_workspace,
      interval_file = interval_file,
      preemptible_tries = preemptible_tries,
      memory_gb = memory_gb,
      disk_gb = disk_gb,
      cpu = cpu,
      gatk_docker = gatk_docker
  }

  output {
    String workspace_dir = GenomicsDBImportTask.workspace_dir
  }
}

task GenomicsDBImportTask {
  meta {
    author: "Shadi Zaheri"
    description: "Uses GATK GenomicsDBImport to consolidate per-sample GVCFs into a GenomicsDB workspace for joint genotyping."
  }

  input {
    File reference_fasta
    File reference_fasta_index
    File reference_dict
    Array[File] gvcf_files
    Array[String] sample_names
    String genomicsdb_workspace
    File? interval_file
    Int preemptible_tries
    Int memory_gb
    Int disk_gb
    Int cpu
    String gatk_docker
  }
  
  command <<<
    set -euo pipefail

    mkdir -p sample_map_dir
    > sample_map_dir/sample_map.txt

  
    for i in $(seq 0 $((${#sample_names[@]} - 1))); do
      echo "${sample_names[$i]}	${gvcf_files[$i]}" >> sample_map_dir/sample_map.txt
    done

    gatk --java-options "-Xmx~{memory_gb}G" GenomicsDBImport \
      --sample-name-map sample_map_dir/sample_map.txt \
      --genomicsdb-workspace-path ~{genomicsdb_workspace} \
      --reference ~{reference_fasta} \
      ~{if defined(interval_file) then "--intervals " + interval_file else ""}
  >>>

  runtime {
    docker: gatk_docker
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: preemptible_tries
    cpu: cpu
  }

  output {
    String workspace_dir = genomicsdb_workspace
  }
}
