version 1.0

workflow BwaAlignmentWorkflow {

    input {
        File fq_end1
        File fq_end2
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        String? read_group
        String prefix = "aligned"

        Int? bwa_cpu_cores
        Int? bwa_memory_gb
        Int? bwa_disk_gb

        Int? sort_cpu_cores
        Int? sort_memory_gb
        Int? sort_disk_gb
    }

    call Bwa {
        input:
            fq_end1 = fq_end1,
            fq_end2 = fq_end2,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            read_group = read_group,
            prefix = prefix,
            cpu_cores = bwa_cpu_cores,
            memory_gb = bwa_memory_gb,
            disk_gb = bwa_disk_gb
    }

    call SortBam {
        input:
            unsorted_bam = Bwa.bam,
            prefix = prefix + ".sorted",
            cpu_cores = sort_cpu_cores,
            memory_gb = sort_memory_gb,
            disk_gb = sort_disk_gb
    }

    output {
        File sorted_bam = SortBam.sorted_bam
        File sorted_bai = SortBam.sorted_bai
    }
}

task Bwa {
    input {
        File fq_end1
        File fq_end2
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        String? read_group    # Optional read group information
        String prefix = "aligned"  # Prefix for output BAM file
        Int? cpu_cores        # Optional number of CPU cores for alignment
        Int? memory_gb        # Optional memory allocation
        Int? disk_gb          # Optional disk space allocation
    }

    Int disk_size = 1 + 4 * ceil(size(fq_end1, "GB")) 
                      + 4 * ceil(size(fq_end2, "GB")) 
                      + 4 * ceil(size(ref_fasta, "GB")) 
                      + 4 * ceil(size(ref_fasta_index, "GB")) 
                      + 4 * ceil(size(ref_dict, "GB"))

    String rg_arg = if defined(read_group) then "-R ~{read_group}" else ""

    command <<<!
        set -euxo pipefail
        # Calculate available processors
        np=$(cat /proc/cpuinfo | grep ^processor | tail -n1 | awk '{print $NF+1}')
        if [[ ${np} -gt 2 ]] ; then
            np=$((np-1))
        fi

        # Perform alignment using bwa mem
        bwa mem \
            -t ${np} \
            ~{rg_arg} \
            ~{ref_fasta} \
            ~{fq_end1} \
            ~{fq_end2} | \
        samtools view -1 - > ~{prefix}.bam
    >>>

    output {
        File bam = "~{prefix}.bam"  # Output unsorted BAM file
    }

    runtime {
        cpu: select_first([cpu_cores, 2])
        memory: select_first([memory_gb, 8]) + " GiB"
        disks: "local-disk " + select_first([disk_gb, disk_size]) + " HDD"
        docker: "us.gcr.io/broad-dsp-lrma/sr-utils:0.2.2"
    }
}

    #########################
    
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/sr-utils:0.2.2"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    runtime {
        cpu:    select_first([runtime_attr.cpu_cores, cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, memory_gb, default_attr.mem_gb]) + " GiB"
        disks:  "local-disk " + select_first([runtime_attr.disk_gb, disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task SortBam {
    input {
        File unsorted_bam
        String prefix = "sorted"
        Int? cpu_cores
        Int? memory_gb
        Int? disk_gb
    }

    Int disk_size = 1 + 2 * ceil(size(unsorted_bam, "GB"))

    command <<<!
        set -euxo pipefail
        samtools sort -@~{cpu_cores} -m ~{memory_gb}G -o ~{prefix}.bam ~{unsorted_bam}
        samtools index ~{prefix}.bam
    >>>

    output {
        File sorted_bam = "~{prefix}.bam"
        File sorted_bai = "~{prefix}.bam.bai"
    }

    runtime {
        cpu: select_first([cpu_cores, 2])
        memory: select_first([memory_gb, 8]) + " GiB"
        disks: "local-disk " + select_first([disk_gb, disk_size]) + " HDD"
        docker: "us.gcr.io/broad-dsp-lrma/sr-utils:0.2.2"
    }
}

