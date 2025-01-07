version 1.0

task SplitVcfs {
    input {
        String donor_id
        String sample_id
        File vcf
        File vcf_index
        String? optional_args = ""
        Int? memory_gb = 16
        Int? disk_size_gb = 100
    }

    Int machine_mem_gb = select_first([memory_gb, 16])
    Int command_mem_gb = machine_mem_gb - 4

    runtime {
        docker: "broadinstitute/gatk:4.2.6.1"
        memory: machine_mem_gb + " GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 3
    }

    command <<<
        set -euo pipefail
        
        gatk --java-options "-Xmx~{command_mem_gb}G -Djava.io.tmpdir=$TMPDIR" \
        SplitVcfs \
        -I ~{vcf} \
        -SNP_OUTPUT ~{sample_id}.snvs.vcf.gz \
        -INDEL_OUTPUT ~{sample_id}.indels.vcf.gz \
        -STRICT false \
        ~{optional_args}
    >>>

    output {
        File snvs_vcf = "~{sample_id}.snvs.vcf.gz"
        File snvs_vcf_index = "~{sample_id}.snvs.vcf.gz.tbi"
        File indels_vcf = "~{sample_id}.indels.vcf.gz"
        File indels_vcf_index = "~{sample_id}.indels.vcf.gz.tbi"
    }
}

workflow SplitVcfsWorkflow {
    input {
        String donor_id
        String sample_id
        File vcf
        File vcf_index
        String? optional_args
    }

    call SplitVcfs {
        input:
            donor_id = donor_id,
            sample_id = sample_id,
            vcf = vcf,
            vcf_index = vcf_index,
            optional_args = optional_args
    }

    output {
        File snvs_vcf = SplitVcfs.snvs_vcf
        File snvs_vcf_index = SplitVcfs.snvs_vcf_index
        File indels_vcf = SplitVcfs.indels_vcf
        File indels_vcf_index = SplitVcfs.indels_vcf_index
    }
}
