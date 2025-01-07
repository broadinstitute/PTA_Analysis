version 1.0

# Import tasks from separate files
import "../../tasks/Utils/createExcludeIndelList.wdl" as Utils
import "../../tasks/bedtools/sort.wdl" as Bedtools

workflow create_excludeindellist_bed {
    input {
        Array[Tuple[String, String, File]] features_beds  # [donor_id, sample_id, feature_bed]
        Array[Tuple[String, String, String]] label_info   # [donor_id, sample_id, label]
    }

    # Filter and prepare input for PTA labels only
    scatter (feature_bed in features_beds) {
        String donor_id = feature_bed.left
        String sample_id = feature_bed.middle
        
        # Find matching label info
        if (donor_id == label_info.left && sample_id == label_info.middle && label_info.right == "PTA") {
            File filtered_beds = feature_bed.right
        }
    }

    # Create exclude indel list
    call Utils.createIndelExcludeIndelList {
        input:
            feature_beds = select_all(filtered_beds)
    }

    # Sort the exclude indel list
    call Bedtools.sort {
        input:
            bed_file = createIndelExcludeIndelList.exclude_indel_list
    }

    # Copy to output directory
    output {
        File excludeindellist_bed = sort.sorted_bed
    }
}