import "../../tasks/randomForest.wdl" as RF

workflow get_indel_rf_files {
    input {
        Array[Tuple[String, String, File, String]] rf_tables_with_labels  # [donor_id, sample_id, rf_table, label]
        String out_dir
    }

    # Group the data by label
    call group_by_label {
        input:
            rf_tables_with_labels = rf_tables_with_labels
    }

    # Train random forest model
    call RF.train_indel_rf {
        input:
            grouped_data = group_by_label.grouped_data
    }

    # Copy files to output directory
    call copy_rf_files {
        input:
            confusion_file = train_indel_rf.confusion_matrix,
            importance_file = train_indel_rf.feature_importance,
            rdata_file = train_indel_rf.rdata,
            rds_file = train_indel_rf.rds_model,
            out_dir = out_dir
    }

    output {
        File output_confusion = copy_rf_files.confusion_out
        File output_importance = copy_rf_files.importance_out
        File output_rdata = copy_rf_files.rdata_out
        File output_rds = copy_rf_files.rds_out
    }
}

task group_by_label {
    input {
        Array[Tuple[String, String, File, String]] rf_tables_with_labels
    }

    command <<<
        # Python script to group data by label
        python <<CODE
        import json
        import sys

        data = ~{write_json(rf_tables_with_labels)}
        grouped = {}
        
        for item in data:
            label = item[3]  # label is the 4th element
            if label not in grouped:
                grouped[label] = []
            grouped[label].append(item[2])  # rf_table is the 3rd element
        
        with open('grouped_data.json', 'w') as f:
            json.dump(grouped, f)
        CODE
    >>>

    output {
        File grouped_data = "grouped_data.json"
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
        mkdir -p ~{out_dir}/indels/randomforest
        cp ~{confusion_file} ~{out_dir}/indels/randomforest/
        cp ~{importance_file} ~{out_dir}/indels/randomforest/
        cp ~{rdata_file} ~{out_dir}/indels/randomforest/
        cp ~{rds_file} ~{out_dir}/indels/randomforest/
    >>>

    output {
        File confusion_out = out_dir + "/indels/randomforest/" + basename(confusion_file)
        File importance_out = out_dir + "/indels/randomforest/" + basename(importance_file)
        File rdata_out = out_dir + "/indels/randomforest/" + basename(rdata_file)
        File rds_out = out_dir + "/indels/randomforest/" + basename(rds_file)
    }

    runtime {
        docker: "ubuntu:latest"
    }
}
