version 1.0

import "../../tasks/randomForest.wdl" as RF

workflow get_snv_rf_files {
    input {
        Array[Tuple[String, String, File, String]] rf_tables_with_labels  # [donor_id, sample_id, rf_table, label]
        String out_dir
    }

    # Group files by label
    call group_by_label {
        input:
            rf_tables_with_labels = rf_tables_with_labels
    }

    # Train random forest model
    call RF.train_snv_rf {
        input:
            grouped_files = group_by_label.grouped_files
    }

    # Copy files to output directory
    call copy_rf_files {
        input:
            confusion_file = train_snv_rf.confusion_matrix,
            importance_file = train_snv_rf.feature_importance,
            rdata_file = train_snv_rf.rdata_file,
            rds_file = train_snv_rf.rds_file,
            out_dir = out_dir
    }

    output {
        File output_confusion = copy_rf_files.final_confusion
        File output_importance = copy_rf_files.final_importance
        File output_rdata = copy_rf_files.final_rdata
        File output_rds = copy_rf_files.final_rds
    }
}

task group_by_label {
    input {
        Array[Tuple[String, String, File, String]] rf_tables_with_labels
    }

    command <<<
        # Python script to group files by label
        python <<CODE
        import json
        import sys

        data = ~{write_json(rf_tables_with_labels)}
        grouped = {}
        
        for item in data:
            label = item[3]  # label is the fourth element
            if label not in grouped:
                grouped[label] = []
            grouped[label].append(item[2])  # rf_table is the third element
        
        with open('grouped_files.json', 'w') as f:
            json.dump(grouped, f)
        CODE
    >>>

    output {
        Map[String, Array[File]] grouped_files = read_json("grouped_files.json")
    }

    runtime {
        docker: "python:3.8-slim"
    }
}

task copy_rf_files {
    input {
        File confusion_file
        File importance_file
        File rdata_file
        File rds_file
        String out_dir
    }

    command <<<
        mkdir -p ~{out_dir}/snvs/randomforest
        cp ~{confusion_file} ~{out_dir}/snvs/randomforest/
        cp ~{importance_file} ~{out_dir}/snvs/randomforest/
        cp ~{rdata_file} ~{out_dir}/snvs/randomforest/
        cp ~{rds_file} ~{out_dir}/snvs/randomforest/
    >>>

    output {
        File final_confusion = out_dir + "/snvs/randomforest/" + basename(confusion_file)
        File final_importance = out_dir + "/snvs/randomforest/" + basename(importance_file)
        File final_rdata = out_dir + "/snvs/randomforest/" + basename(rdata_file)
        File final_rds = out_dir + "/snvs/randomforest/" + basename(rds_file)
    }

    runtime {
        docker: "ubuntu:latest"
    }
}
