version 1.0

import "../../tasks/SMuRF/smurf.wdl" as SMuRF
import "../../tasks/htslib/bgzip.wdl" as Bgzip
import "../../tasks/htslib/tabix.wdl" as Tabix

workflow get_somatic_vcfs {
    input {
        Array[Pair[String, File]] input_germline_vcfs  # [donor_id, vcf]
        Array[Pair[String, File]] input_bams           # [donor_id, bam]
        Array[String] bulk_names
        String out_dir
    }

    # Group inputs by donor_id
    scatter (germline_vcf in input_germline_vcfs) {
        String donor_ids = germline_vcf.left
        
        # Find matching BAM for this donor
        Pair[String, File] matched_bam = select_first([
            input_bams[i] for i in range(length(input_bams)) 
            if input_bams[i].left == donor_ids
        ])
    }

    # Call SMuRF for each donor
    scatter (idx in range(length(donor_ids))) {
        call SMuRF.smurf {
            input:
                donor_id = donor_ids[idx],
                germline_vcf = input_germline_vcfs[idx].right,
                bam = matched_bam[idx].right,
                bulk_name = bulk_names[idx]
        }

        # Compress output VCF
        call Bgzip.bgzip {
            input:
                input_file = smurf.somatic_vcf,
                output_dir = out_dir + "/intermediate/short_variants/somatic_vcfs/" + donor_ids[idx]
        }

        # Index compressed VCF
        call Tabix.tabix {
            input:
                vcf_gz = bgzip.compressed_file
        }
    }

    # Copy intermediate files to output locations
    scatter (idx in range(length(donor_ids))) {
        String sample_id = sub(basename(smurf.somatic_vcf[idx]), ".SMuRF.filtered.vcf$", "")
    }

    output {
        Array[Pair[String, Pair[String, Pair[File, File]]]] somatic_vcfs = zip(
            donor_ids,
            zip(
                sample_id,
                zip(bgzip.compressed_file, tabix.index_file)
            )
        )
    }
}