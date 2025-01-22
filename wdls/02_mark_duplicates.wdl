version 1.0
# Workflow: MarkDuplicatesWorkflow
# Author: Shadi Zaheri
# Date: 2025-01-22
# Description: Marks duplicates in BAM files using GATK MarkDuplicates and provides flexibility for optional parameters.
#
# License: Broad Institute  of MIT and Harvard License
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated 
# documentation files (the "Software"), to deal in the Software without restriction, including without limitation 
# the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, 
# and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or substantial portions 
# of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED 
# TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF 
# CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
# DEALINGS IN THE SOFTWARE.

workflow MarkDuplicatesWorkflow {
   parameter_meta {
    input_bam: "Path to the input BAM file. Must be coordinate-sorted."
    sample_name: "Sample name used for output file naming."
    metrics_filename: "Path to the output duplication metrics file."
    preemptible_tries: "Number of preemptible retries allowed for the task (Google Cloud only)."
    memory_gb: "Memory allocated for the task in gigabytes."
    disk_gb: "Disk space allocated for the task in gigabytes."
    cpu: "Number of CPU cores allocated for the task."
    docker: "Docker image used to run the task (default: broadinstitute/gatk:4.6.1.0)."
    tagging_policy: "Policy for tagging duplicates. Options: All, OpticalOnly, DontTag. Default: All."
    remove_duplicates: "If true, removes duplicate reads entirely instead of marking them. Default: false."
    create_index: "If true, creates an index for the output BAM file. Default: true."
    optical_duplicate_pixel_distance: "Maximum distance in pixels to consider optical duplicates. Default: 100."
    max_file_handles: "Maximum number of file handles to keep open for read ends map. Default: 8000."
    max_optical_duplicate_set_size: "Maximum size of duplicate sets for optical duplicate detection. Default: 300000."
    max_records_in_ram: "Number of records to store in RAM before spilling to disk. Default: 500000."
    sorting_collection_size_ratio: "Fraction of available memory to use for sorting collections. Default: 1.0 (all memory)."
  }

  input {
    File input_bam
    String sample_name
    String metrics_filename
    Int preemptible_tries = 1
    Int memory_gb = 8
    Int disk_gb = 50
    Int cpu = 4
    String docker = "broadinstitute/gatk:4.6.1.0"

    # Optional parameters with user-defined defaults
    String tagging_policy = "All" #"DontTag"    # Options: DontTag, OpticalOnly, All
    Boolean remove_duplicates = false   # Remove duplicates entirely
    Boolean create_index = true         # Create index for the output BAM
    Int optical_duplicate_pixel_distance = 100  # Default for unpatterned flow cells
    Int max_file_handles = 8000         # File handles for read ends map
    Int max_optical_duplicate_set_size = 300000 # Optical duplicate set size
    Int max_records_in_ram = 500000     # Records stored in RAM before spilling to disk
  #  Float sorting_collection_size_ratio = 1.0  # Default to using all available memory
  }

  call MarkDuplicates {
    input:
      input_bam = input_bam,
      sample_name = sample_name,
      metrics_filename = metrics_filename,
      preemptible_tries = preemptible_tries,
      memory_gb = memory_gb,
      disk_gb = disk_gb,
      cpu = cpu,
      docker = docker,
      tagging_policy = tagging_policy,
      remove_duplicates = remove_duplicates,
      create_index = create_index,
      optical_duplicate_pixel_distance = optical_duplicate_pixel_distance,
      max_file_handles = max_file_handles,
      max_optical_duplicate_set_size = max_optical_duplicate_set_size,
      max_records_in_ram = max_records_in_ram
      #sorting_collection_size_ratio = sorting_collection_size_ratio
  }

  output {
    File output_bam = MarkDuplicates.output_bam
    File output_bam_index = MarkDuplicates.output_bam_index
    File duplicate_metrics = MarkDuplicates.duplicate_metrics
  }
}

task MarkDuplicates {
  input {
    File input_bam
    String sample_name
    String metrics_filename
    Int preemptible_tries
    Int memory_gb
    Int disk_gb
    Int cpu
    String docker

    # Optional parameters
    String tagging_policy
    Boolean remove_duplicates
    Boolean create_index
    Int optical_duplicate_pixel_distance
    Int max_file_handles
    Int max_optical_duplicate_set_size
    Int max_records_in_ram
  #  Float sorting_collection_size_ratio
  }

  command {
    set -euxo pipefail

    gatk MarkDuplicates \
      -I ~{input_bam} \
      -O ~{sample_name}.bam \
      -M ~{metrics_filename} \
      --TAGGING_POLICY ~{tagging_policy} \
      --REMOVE_DUPLICATES ~{if remove_duplicates then "true" else "false"} \
      --CREATE_INDEX ~{if create_index then "true" else "false"} \
      --OPTICAL_DUPLICATE_PIXEL_DISTANCE ~{optical_duplicate_pixel_distance} \
      --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP ~{max_file_handles} \
      --MAX_OPTICAL_DUPLICATE_SET_SIZE ~{max_optical_duplicate_set_size} \
      --MAX_RECORDS_IN_RAM ~{max_records_in_ram}
  }

  runtime {
    docker: docker
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_gb} HDD"
    cpu: cpu
    preemptible: preemptible_tries
  }

  output {
    File output_bam = "~{sample_name}.bam"
    File output_bam_index = "~{sample_name}.bam.bai"
    File duplicate_metrics = "~{metrics_filename}"
  }
}
