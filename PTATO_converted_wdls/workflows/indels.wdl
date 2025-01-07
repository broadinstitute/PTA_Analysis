version 1.0

import "indels/get_rf_tables.wdl" as GetRfTables
import "indels/get_excludeindellist_filtered_vcfs.wdl" as GetExcludeIndelListFilteredVcfs
import "indels/get_ptato_vcfs.wdl" as GetPtatoVcfs
import "indels/filter_ptato_vcfs.wdl" as FilterPtatoVcfs
import "indels/get_rf_files.wdl" as GetRfFiles
import "indels/create_excludeindellist_bed.wdl" as CreateExcludeIndelListBed

workflow indels {
    input {
        Array[File] somatic_indel_vcfs
        Array[File] context_beds
    }

    call GetPtatoVcfs.get_ptato_vcfs {
        input:
            somatic_indel_vcfs = somatic_indel_vcfs,
            context_beds = context_beds
    }

    call FilterPtatoVcfs.filter_ptato_vcfs {
        input:
            ptato_vcfs = get_ptato_vcfs.ptato_vcfs
    }

    output {
        Array[File] ptato_vcfs = get_ptato_vcfs.ptato_vcfs
    }
}

workflow indels_train {
    input {
        Array[File] ab_tables
        Array[File] features_beds
        File label_info
    }

    call CreateExcludeIndelListBed.create_excludeindellist_bed {
        input:
            features_beds = features_beds,
            label_info = label_info
    }

    call GetRfTables.get_indel_rf_tables {
        input:
            ab_tables = ab_tables,
            features_beds = features_beds
    }

    call GetRfFiles.get_indel_rf_files {
        input:
            rf_tables = get_indel_rf_tables.rf_tables,
            label_info = label_info
    }
}