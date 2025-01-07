version 1.0

workflow Circos_workflow {
    input {
        String donor_id
        String normal_sample_id
        String tumor_sample_id
        File circos_config_file
        Array[File] circos_txt_files
        String? circos_optional_params = ""  # Optional parameters for circos
    }

    call Circos {
        input:
            donor_id = donor_id,
            normal_sample_id = normal_sample_id,
            tumor_sample_id = tumor_sample_id,
            circos_config_file = circos_config_file,
            circos_txt_files = circos_txt_files,
            circos_optional_params = circos_optional_params
    }

    output {
        Array[File] circos_plots = Circos.circos_plots
    }
}

task Circos {
    input {
        String donor_id
        String normal_sample_id
        String tumor_sample_id
        File circos_config_file
        Array[File] circos_txt_files
        String circos_optional_params

        # Runtime parameters
        String docker = "gridss/gridss-purple-linx:1.3.2"
        Int memory_gb = 8
        Int cpu = 1
        Int disk_size_gb = 50
    }

    command <<<
        set -euo pipefail

        circos \
        -conf ~{circos_config_file} \
        -outputfile ~{tumor_sample_id} \
        ~{circos_optional_params}
    >>>

    output {
        Array[File] circos_plots = glob("${tumor_sample_id}.*g")
    }

    runtime {
        docker: docker
        memory: "~{memory_gb}G"
        cpu: cpu
        disks: "local-disk ~{disk_size_gb} SSD"
    }
}
