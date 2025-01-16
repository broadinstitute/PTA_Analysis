version 1.0

## Copyright Broad Institute, 2018
##
## This WDL defines tasks used for germline variant discovery of human whole-genome or exome sequencing data.
##
## Runtime parameters are often optimized for Broad's Google Cloud Platform implementation.
## For program versions, see docker containers.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

workflow VariantFiltrationWorkflow {

    meta {
        author: "S. Zaheri"
        description: "Filters variants using GATK VariantFiltration with custom filter expressions."
    }

    input {
        File input_vcf       # VCF file from HaplotypeCaller
        File input_vcf_index # VCF index file from HaplotypeCaller
        File ref_fasta       # Reference fasta file
        File ref_fasta_index # Reference fasta index
        File ref_dict        # Reference dictionary
        String sample_name   # Sample name for outputs

        # Optional resource specifications
        Int? cpu             # Number of CPUs
        String? memory       # Memory (e.g., "8 GiB")
        Int? disk_size       # Disk size in GB
    }

    call FilterVariantsGATK {
        input:
            input_vcf = input_vcf,
            input_vcf_index = input_vcf_index,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            sample_name = sample_name,
            runtime_cpu = cpu,
            runtime_memory = memory,
            runtime_disk_size = disk_size
    }

    output {
        File filtered_vcf = FilterVariantsGATK.filtered_vcf
        File filtered_vcf_index = FilterVariantsGATK.filtered_vcf_index
    }
}

task FilterVariantsGATK {

    input {
        File input_vcf       # Input VCF file
        File input_vcf_index # Input VCF index
        File ref_fasta       # Reference fasta file
        File ref_fasta_index # Reference fasta index
        File ref_dict        # Reference dictionary
        String sample_name   # Sample name for outputs

        # Optional runtime inputs
        Int? runtime_cpu
        String? runtime_memory
        Int? runtime_disk_size
    }

    # Default disk size if not provided
    Int disk_size = select_first([runtime_disk_size, ceil(size(input_vcf, "GiB") * 2) + 20])

    command <<<!
        set -euxo pipefail

        gatk --java-options "-Xmx~{select_first([runtime_memory, "4G"])}" \
            VariantFiltration \
            -R ~{ref_fasta} \
            -V ~{input_vcf} \
            --filter-expression "QD < 2.0" --filter-name "SNP_LowQualityDepth" \
            --filter-expression "MQ < 40.0" --filter-name "SNP_MappingQuality" \
            --filter-expression "FS > 60.0" --filter-name "SNP_StrandBias" \
            --filter-expression "HaplotypeScore > 13.0" --filter-name "SNP_HaplotypeScoreHigh" \
            --filter-expression "MQRankSum < -12.5" --filter-name "SNP_MQRankSumLow" \
            --filter-expression "ReadPosRankSum < -8.0" --filter-name "SNP_ReadPosRankSumLow" \
            --filter-expression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" --filter-name "SNP_HardToValidate" \
            --filter-expression "DP < 5" --filter-name "SNP_LowCoverage" \
            --filter-expression "QUAL < 30" --filter-name "SNP_VeryLowQual" \
            --filter-expression "QUAL >= 30.0 && QUAL < 50.0" --filter-name "SNP_LowQual" \
            --filter-expression "SOR > 4.0" --filter-name "SNP_SOR" \
            --cluster 3 \
            --window 10 \
            -O ~{sample_name}.filtered.vcf.gz
    >>>

    runtime {
        docker: "us.gcr.io/broad-gatk/gatk:latest"
        memory: select_first([runtime_memory, "6 GiB"])
        cpu: select_first([runtime_cpu, 2])
        disks: "local-disk " + disk_size + " HDD"
    }

    output {
        File filtered_vcf = "~{sample_name}.filtered.vcf.gz"
        File filtered_vcf_index = "~{sample_name}.filtered.vcf.gz.tbi"
    }
}
