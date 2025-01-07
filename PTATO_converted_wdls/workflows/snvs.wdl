version 1.0

import "./snvs/get_rf_tables.wdl" as get_rf_tables
import "./snvs/get_ptato_vcfs.wdl" as get_ptato_vcfs
import "./snvs/filter_ptato_vcfs.wdl" as filter_ptato_vcfs
import "./snvs/get_rf_files.wdl" as get_rf_files

workflow snvs {
    input {
        Array[File] ab_tables
        Array[File] features_beds
        Array[File] somatic_snv_vcfs
        Array[File] walker_vcfs
    }

    call get_rf_tables.get_snvs_rf_tables {
        input:
            ab_tables = ab_tables,
            features_beds = features_beds
    }

    call get_ptato_vcfs.get_ptato_vcfs {
        input:
            somatic_snv_vcfs = somatic_snv_vcfs,
            rf_tables = get_snvs_rf_tables.rf_tables
    }

    call filter_ptato_vcfs.filter_ptato_vcfs {
        input:
            ptato_vcfs = get_ptato_vcfs.ptato_vcfs,
            walker_vcfs = walker_vcfs
    }

    output {
        Array[Pair[File, File]] ptato_combined_vcfs = zip(get_ptato_vcfs.ptato_vcfs, filter_ptato_vcfs.ptato_filtered_vcfs)
    }
}

workflow snvs_train {
    input {
        Array[File] ab_tables
        Array[File] features_beds
        File label_info
    }

    call get_rf_tables.get_snvs_rf_tables {
        input:
            ab_tables = ab_tables,
            features_beds = features_beds
    }

    call get_rf_files.get_snv_rf_files {
        input:
            rf_tables = get_snvs_rf_tables.rf_tables,
            label_info = label_info
    }

    output {
        Array[File] rf_files = get_snv_rf_files.rf_files
    }
}
