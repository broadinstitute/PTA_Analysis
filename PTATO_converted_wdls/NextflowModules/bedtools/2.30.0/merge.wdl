version 1.0

workflow bedtools_merge {
    input {
        # Workflow level inputs
        String donor_id
        String sample_id
        Array[File] feature_beds
        Array[String] feature_ids
        String? merge_params = ""  # Optional parameter with empty default
    }

    # Scatter over feature beds to merge each feature separately
    scatter (pair in zip(feature_beds, feature_ids)) {
        call merge {
            input:
                donor_id = donor_id,
                sample_id = sample_id,
                feature_id = pair.right,
                bed = pair.left,
                merge_params = merge_params
        }
    }

    # Merge all features together
    call mergeAll {
        input:
            donor_id = donor_id,
            sample_id = sample_id,
            bed = feature_beds[0]  # You might want to merge all feature_beds first
    }

    output {
        Array[File] feature_merged_beds = merge.feature_merged_bed
        File all_features_merged = mergeAll.feature_merged_bed
    }
}

task merge {
    input {
        String donor_id
        String sample_id
        String feature_id
        File bed
        String merge_params
        # You might want to make these runtime parameters configurable
        Int memory_gb = 4
        Int disk_size_gb = 10
        Int cpu = 1
    }

    command <<<
        hostname
        
        bedtools \
        merge \
        -i ~{bed} \
        ~{merge_params} \
        > ~{sample_id}.~{feature_id}.merged.bed
    >>>

    output {
        File feature_merged_bed = "~{sample_id}.~{feature_id}.merged.bed"
    }

    runtime {
        docker: "quay.io/biocontainers/bedtools:2.30.0--h468198e_3"
        memory: "~{memory_gb}G"
        disks: "local-disk ~{disk_size_gb} SSD"
        cpu: cpu
    }
}

task mergeAll {
    input {
        String donor_id
        String sample_id
        File bed
        # You might want to make these runtime parameters configurable
        Int memory_gb = 4
        Int disk_size_gb = 10
        Int cpu = 1
    }

    command <<<
        hostname
        
        bedtools \
        merge \
        -i ~{bed} \
        > ~{sample_id}.features.merged.bed
    >>>

    output {
        File feature_merged_bed = "~{sample_id}.features.merged.bed"
    }

    runtime {
        docker: "quay.io/biocontainers/bedtools:2.30.0--h468198e_3"
        memory: "~{memory_gb}G"
        disks: "local-disk ~{disk_size_gb} SSD"
        cpu: cpu
    }
}
