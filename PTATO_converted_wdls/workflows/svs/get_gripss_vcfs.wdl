version 1.0

# Import the tasks that correspond to GripssApplicationKt and GripssHardFilterApplicationKt
import "../../tasks/gripss/gripss.wdl" as gripss

workflow get_gripss_vcfs {
    input {
        # Input is a tuple of [donor_id, normal_sample_id, tumor_sample_id, gridss_vcf, gridss_tbi]
        Array[Tuple[String, String, String, File, File]] gridss_unfiltered_vcfs
        String output_dir
    }

    scatter (vcf_tuple in gridss_unfiltered_vcfs) {
        call gripss.run_gripss {
            input:
                donor_id = vcf_tuple.0,
                normal_sample_id = vcf_tuple.1,
                tumor_sample_id = vcf_tuple.2,
                input_vcf = vcf_tuple.3,
                input_vcf_index = vcf_tuple.4,
                output_dir = output_dir + "/intermediate/svs/gripss/" + vcf_tuple.0 + "/" + vcf_tuple.1
        }

        call gripss.run_gripss_hard_filter {
            input:
                donor_id = vcf_tuple.0,
                normal_sample_id = vcf_tuple.1,
                tumor_sample_id = vcf_tuple.2,
                gripss_vcf = run_gripss.output_vcf,
                gripss_vcf_index = run_gripss.output_vcf_index,
                output_dir = output_dir + "/intermediate/svs/gripss/" + vcf_tuple.0 + "/" + vcf_tuple.1
        }
    }

    output {
        Array[Tuple[String, String, String, File, File]] gripss_somatic_filtered_vcfs = zip(
            gridss_unfiltered_vcfs.0,  # donor_id
            gridss_unfiltered_vcfs.1,  # normal_sample_id
            gridss_unfiltered_vcfs.2,  # tumor_sample_id
            run_gripss_hard_filter.output_vcf,
            run_gripss_hard_filter.output_vcf_index
        )
    }
}
