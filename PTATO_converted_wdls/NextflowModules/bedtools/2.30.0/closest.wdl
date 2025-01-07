version 1.0

workflow bedtools_closest {
    input {
        String donor_id
        String sample_id
        File bed
        String feature_id
        File feature_bed
        String merge_params
        String? optional_params
    }

    call closest {
        input:
            donor_id = donor_id,
            sample_id = sample_id,
            bed = bed,
            feature_id = feature_id,
            feature_bed = feature_bed,
            merge_params = merge_params,
            optional_params = optional_params
    }

    output {
        File output_bed = closest.feature_bed
    }
}

task closest {
  input {
    String donor_id
    String sample_id
    File bed
    String feature_id
    File feature_bed
    String merge_params
    String? optional_params # Equivalent to params.bedtoolsclosest.optional
  }

  command <<<
    hostname
    
    bedtools \
    closest \
    -a ~{bed} \
    -b ~{feature_bed} \
    ~{optional_params} \
    > ~{sample_id}.~{feature_id}.bed
  >>>

  output {
    File feature_bed = "~{sample_id}.~{feature_id}.bed"
  }

  runtime {
    docker: "quay.io/biocontainers/bedtools:2.30.0--h468198e_3"
  }
}
