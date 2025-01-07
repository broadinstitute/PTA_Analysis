version 1.0

# Import tasks (these would need to be defined in separate files)
import "../../tasks/getContext.wdl" as GetContextTask
import "../../tasks/bedtools.wdl" as BedtoolsTask

workflow get_context_beds {
    input {
        Array[Tuple[String, String, File]] vcfs  # Array of [donor_id, sample_id, vcf_file]
    }

    scatter (vcf_tuple in vcfs) {
        String donor_id = vcf_tuple.left
        String sample_id = vcf_tuple.middle
        File vcf_file = vcf_tuple.right

        call GetContextTask.getContext {
            input:
                vcf = vcf_file
        }

        call BedtoolsTask.sort {
            input:
                bed = getContext.bed_out
        }

        # Rename and move the file to output directory
        String bed_name = sub(basename(sort.sorted_bed), ".sorted.bed", ".context.sorted.bed")
        String output_path = "intermediate/short_variants/context/" + donor_id + "/" + bed_name
    }

    output {
        Array[Tuple[String, String, File]] context_beds = zip(
            donor_id,
            sample_id,
            sort.sorted_bed
        )
    }
}