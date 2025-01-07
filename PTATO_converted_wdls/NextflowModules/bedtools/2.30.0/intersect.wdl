version 1.0

workflow bedtools_intersect {
    input {
        String donor_id
        String sample_id
        
        # Basic intersect inputs
        File? input_bed
        Array[File] feature_beds = []
        Array[String] feature_ids = []
        String? merge_params = ""
        
        # IntersectAll inputs
        Array[File] closest_feature_beds = []
        Array[File] intersect_feature_beds = []
        
        # PTATO-specific inputs
        File? input_vcf
        File? input_tbi
        Array[String] ptato_snvs_sample_ids = []
        Array[File] ptato_snvs_vcfs = []
        Array[File] ptato_snvs_tbis = []
        Array[String] ptato_indels_sample_ids = []
        Array[File] ptato_indels_vcfs = []
        Array[File] ptato_indels_tbis = []
        
        # Indel exclusion inputs
        File? excludeindellist
    }

    # Basic intersect operations for each feature
    if (defined(input_bed) && length(feature_beds) > 0) {
        scatter (idx in range(length(feature_beds))) {
            call intersect {
                input:
                    donor_id = donor_id,
                    sample_id = sample_id,
                    bed = select_first([input_bed]),
                    feature_id = feature_ids[idx],
                    feature_bed = feature_beds[idx],
                    merge_params = merge_params
            }
        }
    }

    # Intersect all features together
    if (defined(input_bed) && (length(closest_feature_beds) > 0 || length(intersect_feature_beds) > 0)) {
        call intersectAll {
            input:
                donor_id = donor_id,
                sample_id = sample_id,
                bed = select_first([input_bed]),
                closest_feature_beds = closest_feature_beds,
                intersect_feature_beds = intersect_feature_beds
        }
    }

    # PTATO-specific intersections
    if (defined(input_vcf) && defined(input_tbi) && 
        (length(ptato_snvs_vcfs) > 0 || length(ptato_indels_vcfs) > 0)) {
        call intersectPTATO {
            input:
                donor_id = donor_id,
                sample_id = sample_id,
                input_vcf = select_first([input_vcf]),
                input_tbi = select_first([input_tbi]),
                ptato_snvs_sample_ids = ptato_snvs_sample_ids,
                ptato_snvs_vcfs = ptato_snvs_vcfs,
                ptato_snvs_tbis = ptato_snvs_tbis,
                ptato_indels_sample_ids = ptato_indels_sample_ids,
                ptato_indels_vcfs = ptato_indels_vcfs,
                ptato_indels_tbis = ptato_indels_tbis
        }
    }

    # Indel exclusion list filtering
    if (defined(input_vcf) && defined(input_tbi) && defined(excludeindellist)) {
        call intersectExcludeIndelList {
            input:
                donor_id = donor_id,
                sample_id = sample_id,
                vcf = select_first([input_vcf]),
                tbi = select_first([input_tbi]),
                excludeindellist = select_first([excludeindellist])
        }
    }

    output {
        Array[File]? feature_beds = intersect.feature_bed
        File? features_bed = intersectAll.features_bed
        File? ptato_intersect_vcf = intersectPTATO.ptato_intersect_vcf
        File? excludeindellist_filtered_vcf = intersectExcludeIndelList.excludeindellist_filtered_vcf
    }
}

task intersect {
    input {
        String donor_id
        String sample_id
        File bed
        String feature_id
        File feature_bed
        String merge_params
        # You'll need to define these params in your workflow
        String bedtoolsintersect_optional = ""
        String docker = "quay.io/biocontainers/bedtools:2.30.0--h468198e_3"
    }

    command <<<
        hostname
        
        bedtools \
        intersect \
        -a ~{bed} \
        -b ~{feature_bed} \
        ~{bedtoolsintersect_optional} \
        > ~{sample_id}.~{feature_id}.bed
    >>>

    output {
        File feature_bed = "~{sample_id}.~{feature_id}.bed"
    }

    runtime {
        docker: docker
        shell: "/bin/bash -euo pipefail"
    }
}

task intersectAll {
    input {
        String donor_id
        String sample_id
        File bed
        Array[File] closest_feature_beds = []
        Array[File] intersect_feature_beds = []
        String bedtoolsintersectall_optional = ""
        String docker = "quay.io/biocontainers/bedtools:2.30.0--h468198e_3"
    }

    command <<<
        hostname
        names=""
        regex=".+\.([A-Z]+).groupby.bed"
        
        for feature in ~{sep=' ' closest_feature_beds} ~{sep=' ' intersect_feature_beds}; do
            if [[ $feature =~ $regex ]]; then
                name="${BASH_REMATCH[1]}"
                names=$names' '$name
            fi
        done

        bedtools \
        intersect \
        -a ~{bed} \
        ~{true='-b ' false='' length(closest_feature_beds) > 0}~{sep=' -b ' closest_feature_beds} \
        ~{true='-b ' false='' length(intersect_feature_beds) > 0}~{sep=' -b ' intersect_feature_beds} \
        -names ${names} \
        ~{bedtoolsintersectall_optional} \
        > ~{sample_id}.features.bed
    >>>

    output {
        File features_bed = "~{sample_id}.features.bed"
    }

    runtime {
        docker: docker
        shell: "/bin/bash -euo pipefail"
    }
}

task intersectExcludeIndelList {
    input {
        String donor_id
        String sample_id
        File vcf
        File tbi
        File excludeindellist
        String bedtoolsintersectexcludeindellist_optional = ""
        String docker = "quay.io/biocontainers/bedtools:2.30.0--h468198e_3"
    }

    command <<<
        hostname

        bedtools \
        intersect \
        -a ~{vcf} \
        -b ~{excludeindellist} \
        ~{bedtoolsintersectexcludeindellist_optional} \
        > ~{sample_id}.excludeindellist.filtered.vcf
    >>>

    output {
        File excludeindellist_filtered_vcf = "~{sample_id}.excludeindellist.filtered.vcf"
    }

    runtime {
        docker: docker
        shell: "/bin/bash -euo pipefail"
    }
}

task intersectPTATO {
    input {
        String donor_id
        String sample_id
        File input_vcf
        File input_tbi
        Array[String] ptato_snvs_sample_ids = []
        Array[File] ptato_snvs_vcfs = []
        Array[File] ptato_snvs_tbis = []
        Array[String] ptato_indels_sample_ids = []
        Array[File] ptato_indels_vcfs = []
        Array[File] ptato_indels_tbis = []
        String bedtoolsintersectptato_optional = ""
        String docker = "quay.io/biocontainers/bedtools:2.30.0--h468198e_3"
    }

    command <<<
        hostname

        bedtools \
        intersect \
        -a ~{input_vcf} \
        ~{true='-b ' false='' length(ptato_snvs_vcfs) > 0}~{sep=' -b ' ptato_snvs_vcfs} \
        ~{true='-b ' false='' length(ptato_indels_vcfs) > 0}~{sep=' -b ' ptato_indels_vcfs} \
        ~{bedtoolsintersectptato_optional} \
        > ~{donor_id}.ptato.intersect.vcf
    >>>

    output {
        File ptato_intersect_vcf = "~{donor_id}.ptato.intersect.vcf"
    }

    runtime {
        docker: docker
        shell: "/bin/bash -euo pipefail"
    }
}