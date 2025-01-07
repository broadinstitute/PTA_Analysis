version 1.0

task GetSampleName {
    input {
        String donor_id
        String sample_id
        File bam
        File bai
        String? optional_args
        Int memory_gb = 16
        Int disk_size_gb = ceil(size(bam, "GB") * 2 + 20)
    }

    command <<<
        gatk --java-options "-Xmx~{memory_gb-4}G" \
        GetSampleName \
        -I ~{bam} \
        -O /dev/stdout \
        ~{optional_args}
    >>>

    runtime {
        docker: "broadinstitute/gatk:4.2.6.1"
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_size_gb} HDD"
        preemptible: 3
    }

    output {
        String sample_name = read_string(stdout())
    }
}

workflow GetSampleName_wf {
    input {
        String donor_id
        String sample_id
        File bam
        File bai
        String? optional_args
    }

    call GetSampleName {
        input:
            donor_id = donor_id,
            sample_id = sample_id,
            bam = bam,
            bai = bai,
            optional_args = optional_args
    }

    output {
        String sample_name = GetSampleName.sample_name
    }
}
