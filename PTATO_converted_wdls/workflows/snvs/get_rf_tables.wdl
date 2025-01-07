version 1.0

import "../../tasks/Utils/getFilesFromDir.wdl" as Utils
import "../../tasks/Utils/createRfTable.wdl" as RfTable

workflow get_snvs_rf_tables {
    input {
        Array[Pair[String, Pair[String, File]]] ab_tables    # [donor_id, [sample_id, file]]
        Array[Pair[String, Pair[String, File]]] features_beds # [donor_id, [sample_id, file]]
        String? rf_tables_dir
        String out_dir
    }

    if (defined(rf_tables_dir)) {
        call Utils.extractSnvRfTableRdsFromDir {
            input:
                directory = select_first([rf_tables_dir])
        }
    }

    if (!defined(rf_tables_dir)) {
        scatter (idx in range(length(ab_tables))) {
            String donor_id = ab_tables[idx].left
            String sample_id = ab_tables[idx].right.left
            File ab_table = ab_tables[idx].right.right
            File features_bed = features_beds[idx].right.right

            call RfTable.createSnvRfTable {
                input:
                    donor_id = donor_id,
                    sample_id = sample_id,
                    ab_table = ab_table,
                    features_bed = features_bed,
                    out_dir = out_dir + "/intermediate/snvs/rf/" + donor_id
            }
        }
    }

    output {
        Array[Pair[String, Pair[String, File]]] rf_tables = select_first([
            extractSnvRfTableRdsFromDir.rf_tables,
            zip(
                flatten([donor_id]),
                zip(
                    flatten([sample_id]),
                    createSnvRfTable.rf_table_rds
                )
            )
        ])
    }
}
