version 1.0

task bwa_mem2_alignment {
    input {
        String fastq1
        String fastq2
        Array[File] bwaIndexFiles  # BWA index files
        String reference
        String prefix
        Int threads
        String read_group
        String memory
        String disk
    }

    command {
        echo "Running bwa-mem2..."
        bwa-mem2 mem \
            -K 100000000 \
            -v 3 \
            -t ~{threads} \
            -Y \
            -R "~{read_group}" \
            -M \
            ~{reference} \
            ~{fastq1} ~{fastq2} | \
            samtools view -b -o ~{prefix}.bam -
    }

    output {
        File bam = "~{prefix}_unsorted.bam"
    }

    runtime {
        cpu: threads
        memory: memory
        disks: "local-disk ~{disk}" #HDD or SSD
        docker: "us.gcr.io/broad-dsp-lrma/sr-utils:0.2.2"
    }
}

workflow alignment_workflow {
    input {
        String fastq1
        String fastq2
        Array[File] bwaIndexFiles  # Input the BWA index files
        String reference
        String prefix
        Int threads = 16
        String read_group
        String memory = "16 GB"
        String disk = "200 GB"
    }

    call bwa_mem2_alignment {
        input:
            fastq1 = fastq1,
            fastq2 = fastq2,
            bwaIndexFiles = bwaIndexFiles,  # Pass the index files to the task
            reference = reference,
            prefix = prefix,
            threads = threads,
            read_group = read_group,
            memory = memory,
            disk = disk
    }

    output {
        File bam = bwa_mem2_alignment.bam
    }
}
