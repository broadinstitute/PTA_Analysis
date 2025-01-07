version 1.0

workflow ABtable {
    input {
        Array[Tuple[String, String, String, File, File, File, File, File, File]] input_tuples
        # Each tuple contains: donor_id, sample_id, chrom, phased_vcf, phased_tbi, germline_vcf, germline_tbi, somatic_vcf, somatic_tbi
    }

    scatter (input_tuple in input_tuples) {
        call createABtable {
            input:
                donor_id = input_tuple.0,
                sample_id = input_tuple.1,
                chrom = input_tuple.2,
                phased_vcf = input_tuple.3,
                phased_tbi = input_tuple.4,
                germline_vcf = input_tuple.5,
                germline_tbi = input_tuple.6,
                somatic_vcf = input_tuple.7,
                somatic_tbi = input_tuple.8
        }
    }

    call mergeABtable {
        input:
            donor_id = input_tuples[0].0,  # Taking first donor_id as they should all be the same
            sample_id = input_tuples[0].1,  # Taking first sample_id as they should all be the same
            ab_tables = createABtable.ab_table_out
    }

    output {
        File final_ab_table = mergeABtable.ab_table_out
    }
}

task createABtable {
    input {
        String donor_id
        String sample_id
        String chrom
        File phased_vcf
        File phased_tbi
        File germline_vcf
        File germline_tbi
        File somatic_vcf
        File somatic_tbi
    }

    command <<<
        hostname > host.txt
        Rscript ~{baseDir}/scripts/R/ABscript.R \
            ~{somatic_vcf} \
            ~{germline_vcf} \
            ~{phased_vcf} \
            ~{chrom} \
            ~{sample_id}_~{chrom}.abtable.txt \
            ~{ref_genome}
    >>>

    output {
        File ab_table_out = "~{sample_id}_~{chrom}.abtable.txt"
    }

    runtime {
        docker: "your_docker_image"
        memory: "4 GB"
        cpu: 1
    }
}

task mergeABtable {
    input {
        String donor_id
        String sample_id
        Array[File] ab_tables
    }

    command <<<
        header=true
        for ab_table in ~{sep=' ' ab_tables}; do
            if $header; then
                cat $ab_table > ~{sample_id}.abtable.txt
                header=false
            else
                tail -n+2 $ab_table >> ~{sample_id}.abtable.txt
            fi
        done
    >>>

    output {
        File ab_table_out = "~{sample_id}.abtable.txt"
    }

    runtime {
        docker: "your_docker_image"
        memory: "4 GB"
        cpu: 1
    }
}
