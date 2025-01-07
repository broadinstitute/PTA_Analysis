task CallableLoci {
    input {
        String donor_id
        String sample_id
        File bam
        File bai
        File ref_fasta
        String? optional_args
        Int memory_gb = 16
        Int disk_size_gb = ceil(size(bam, "GB") * 2 + size(ref_fasta, "GB")) + 20
    }

    command <<<
        java -Djava.io.tmpdir=${tmp} \
        -Xmx~{memory_gb-4}G -jar \
        /usr/local/opt/gatk-3.8/GenomeAnalysisTK.jar \
        -T CallableLoci \
        -I ~{bam} \
        -R ~{ref_fasta} \
        -o ~{sample_id}.callableloci.bed \
        --summary ~{sample_id}.callableloci.txt \
        ~{optional_args}
    >>>

    runtime {
        docker: "quay.io/biocontainers/gatk:3.8--hdfd78af_11"
        memory: "~{memory_gb} GB"
        cpu: 1
        disks: "local-disk ~{disk_size_gb} SSD"
    }

    output {
        File callable_loci_bed = "~{sample_id}.callableloci.bed"
        File callable_loci_summary = "~{sample_id}.callableloci.txt"
    }

    meta {
        author: "Converted from Nextflow"
        description: "Runs GATK CallableLoci to determine callable regions in a BAM file"
    }
}

workflow GATK_CallableLoci {
    input {
        String donor_id
        String sample_id
        File bam
        File bai
        File ref_fasta
        String? optional_args
    }

    call CallableLoci {
        input:
            donor_id = donor_id,
            sample_id = sample_id,
            bam = bam,
            bai = bai,
            ref_fasta = ref_fasta,
            optional_args = optional_args
    }

    output {
        File callable_loci_bed = CallableLoci.callable_loci_bed
        File callable_loci_summary = CallableLoci.callable_loci_summary
    }
}
