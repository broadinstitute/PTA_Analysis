version 1.0

task groupby {
    input {
        String donor_id
        String sample_id
        String feature_id
        File bed_file
        String groupby_params
        String? optional_params = ""  # equivalent to params.bedtoolsgroupby.optional
    }

    command <<<
        hostname
        
        bedtools \
        groupby \
        -i ~{bed_file} \
        ~{groupby_params} \
        ~{optional_params} \
        > ~{sample_id}.~{feature_id}.groupby.bed
    >>>

    output {
        File feature_groupby_bed = "~{sample_id}.~{feature_id}.groupby.bed"
    }

    runtime {
        docker: "quay.io/biocontainers/bedtools:2.30.0--h468198e_3"
    }
}

task groupbyAll {
    input {
        String donor_id
        String sample_id
        File bed_file
        String? optional_params = ""  # equivalent to params.bedtoolsgroupbyall.optional
    }

    command <<<
        hostname
        
        bedtools \
        groupby \
        -i ~{bed_file} \
        ~{optional_params} \
        > ~{sample_id}.features.groupby.bed
    >>>

    output {
        File feature_groupby_bed = "~{sample_id}.features.groupby.bed"
    }

    runtime {
        docker: "quay.io/biocontainers/bedtools:2.30.0--h468198e_3"
    }
}

workflow bedtools_groupby {
    input {
        # Inputs for groupby task
        String donor_id
        String sample_id
        String feature_id
        File bed_file
        String groupby_params
        
        # Whether to run groupbyAll instead of groupby
        Boolean run_groupby_all = false
    }

    if (!run_groupby_all) {
        call groupby {
            input:
                donor_id = donor_id,
                sample_id = sample_id,
                feature_id = feature_id,
                bed_file = bed_file,
                groupby_params = groupby_params
        }
    }

    if (run_groupby_all) {
        call groupbyAll {
            input:
                donor_id = donor_id,
                sample_id = sample_id,
                bed_file = bed_file
        }
    }

    output {
        File? groupby_output = groupby.feature_groupby_bed
        File? groupby_all_output = groupbyAll.feature_groupby_bed
    }
}