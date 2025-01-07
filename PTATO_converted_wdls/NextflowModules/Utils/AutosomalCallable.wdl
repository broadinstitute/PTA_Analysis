version 1.0

workflow AutosomalCallableWorkflow {
    input {
        String donor_id
        String sample_id
        File callable_bed_file
        File callable_txt_file
    }

    call AutosomalCallableLoci {
        input:
            donor_id = donor_id,
            sample_id = sample_id,
            callable_bed_file = callable_bed_file,
            callable_txt_file = callable_txt_file
    }

    output {
        File autosomal_callable_file = AutosomalCallableLoci.autosomal_callable_file
    }
}

task AutosomalCallableLoci {
    input {
        String donor_id
        String sample_id
        File callable_bed_file
        File callable_txt_file
    }

    command <<<
        set -euo pipefail
        HOSTNAME=$(hostname)
        echo "Running on host: $HOSTNAME"
        
        Rscript ~{default="/scripts/R/AutosomalCallableLoci.R"} \
            ~{callable_bed_file} \
            ~{sample_id}.callableloci.autosomal.txt
    >>>

    output {
        File autosomal_callable_file = "~{sample_id}.callableloci.autosomal.txt"
        String hostname = read_string(stdout())  # If you want to capture the hostname as output
    }

    runtime {
        docker: "rocker/r-base:latest"  # Example R docker image
        memory: "4 GB"
        cpu: 1
        disks: "local-disk 10 SSD"
    }
}