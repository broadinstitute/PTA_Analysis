version 1.0

workflow filter_ptato_vcfs {
    input {
        Array[Pair[String, Pair[String, File]]] ptato_vcfs    # [donor_id, [sample_id, vcf]]
        Array[Pair[String, Pair[String, File]]] walker_vcfs   # [donor_id, [sample_id, vcf]]
        String out_dir
    }

    # Match up corresponding VCFs by donor_id and sample_id
    scatter (ptato_vcf in ptato_vcfs) {
        String donor_id = ptato_vcf.left
        String sample_id = ptato_vcf.right.left
        File ptato_file = ptato_vcf.right.right

        # Find matching walker VCF
        File walker_file = select_first([
            walker_vcfs[x].right.right 
            if walker_vcfs[x].left == donor_id && walker_vcfs[x].right.left == sample_id
        ])

        call ptato_cutoff {
            input:
                donor_id = donor_id,
                sample_id = sample_id,
                ptato_vcf = ptato_file,
                walker_vcf = walker_file
        }

        call ptato_filter {
            input:
                donor_id = donor_id,
                sample_id = sample_id,
                ptato_vcf = ptato_file,
                walker_vcf = walker_file,
                ptaprob_cutoff = ptato_cutoff.cutoff_value
        }

        call bgzip {
            input:
                vcf = ptato_filter.filtered_vcf,
                output_filename = "${donor_id}.${sample_id}.filtered.vcf.gz"
        }

        call tabix {
            input:
                vcf_gz = bgzip.zipped_vcf
        }
    }

    output {
        Array[Tuple[String, String, File, File, File]] filtered_results = zip(
            donor_id,
            sample_id,
            bgzip.zipped_vcf,
            tabix.vcf_index,
            ptato_cutoff.ptato_table
        )
    }
}

task ptato_cutoff {
    input {
        String donor_id
        String sample_id
        File ptato_vcf
        File walker_vcf
    }

    command <<<
        # Insert actual ptato_cutoff command here
        ptato_cutoff ~{ptato_vcf} ~{walker_vcf}
    >>>

    output {
        File ptato_table = "ptato_table.txt"
        Float cutoff_value = read_float("cutoff.txt")
    }
}

task ptato_filter {
    input {
        String donor_id
        String sample_id
        File ptato_vcf
        File walker_vcf
        Float ptaprob_cutoff
    }

    command <<<
        # Insert actual ptato_filter command here
        ptato_filter ~{ptato_vcf} ~{walker_vcf} ~{ptaprob_cutoff}
    >>>

    output {
        File filtered_vcf = "filtered.vcf"
    }
}

task bgzip {
    input {
        File vcf
        String output_filename
    }

    command <<<
        bgzip -c ~{vcf} > ~{output_filename}
    >>>

    output {
        File zipped_vcf = output_filename
    }
}

task tabix {
    input {
        File vcf_gz
    }

    command <<<
        tabix -p vcf ~{vcf_gz}
    >>>

    output {
        File vcf_index = "~{vcf_gz}.tbi"
    }
}