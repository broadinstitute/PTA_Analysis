version 1.0

import "../../tasks/Utils/getFilesFromDir.wdl" as Utils
import "../../tasks/Utils/createRfTable.wdl" as RfTable

workflow get_indel_rf_tables {
    input {
        Array[Tuple[String, String, File]] ab_tables
        Array[Tuple[String, String, File]] features_beds
        String? rf_tables_dir
        String out_dir
    }

    if (defined(rf_tables_dir)) {
        call Utils.extractIndelRfTableRdsFromDir {
            input:
                dir = select_first([rf_tables_dir])
        }
    }

    if (!defined(rf_tables_dir)) {
        scatter (pair in cross(ab_tables, features_beds)) {
            if (pair.left[0] == pair.right[0] && pair.left[1] == pair.right[1]) {
                call RfTable.createIndelRfTable {
                    input:
                        ab_table = pair.left[2],
                        features_bed = pair.right[2],
                        donor_id = pair.left[0],
                        sample_id = pair.left[1],
                        out_dir = out_dir + "/intermediate/indels/rf/" + pair.left[0]
                }
            }
        }
    }

    output {
        Array[Tuple[String, String, File]] rf_tables = select_first([
            extractIndelRfTableRdsFromDir.rf_tables,
            select_all(createIndelRfTable.rf_table_out)
        ])
    }
}
