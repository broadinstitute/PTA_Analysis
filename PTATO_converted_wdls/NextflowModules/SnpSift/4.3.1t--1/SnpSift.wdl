workflow SnpSiftWorkflow {
  input {
    String donor_id
    String sample_id
    File vcf_gz
    File tbi
    Array[String] normal_sample_ids
    String? snpsift_optional = ""  # Optional parameter with a default empty string
  }

  call SnpSift {
    input:
      donor_id = donor_id,
      sample_id = sample_id,
      vcf_gz = vcf_gz,
      tbi = tbi,
      normal_sample_ids = normal_sample_ids,
      params.snpsift.optional = snpsift_optional
  }

  output {
    File phased_vcf = SnpSift.phased_vcf
  }
}

task SnpSift {
  input {
    String donor_id
    String sample_id
    File vcf_gz
    File tbi
    Array[String] normal_sample_ids
    String snpsift_optional
  }

  command <<<
    host=$(hostname)
    echo ${host}

    isHET_ARRAY=()
    VAF_ARRAY=()
    HEADER=$(zcat ${vcf_gz} | grep "^#CHROM" | cut -f 10-)
    N=$(echo ${HEADER} | grep -o " " | wc -l)

    for NORMAL in ${normal_sample_ids[@]}; do
      NORMAL_GEN_ID=$(echo ${HEADER/${NORMAL}//} | cut -d/ -f1 | wc -w | tr -d ' ')
      isHET="isHet(GEN[${NORMAL_GEN_ID}])"
      VAF="(GEN[${NORMAL_GEN_ID}].AD[1]/(GEN[${NORMAL_GEN_ID}].AD[0]+0.001+GEN[${NORMAL_GEN_ID}].AD[1]))>0.3"
      isHET_ARRAY+=(${isHET})
      VAF_ARRAY+=(${VAF})
    done
    isHET_STRING=${isHET_ARRAY[@]}
    isHet=${isHET_STRING// / & }
    VAF_STRING=${VAF_ARRAY[@]}
    vaf=${VAF_STRING// / & }

    zcat ${vcf_gz} | \
    SnpSift filter \
    "( (exists ID) & (ID =~ 'rs') & ${isHet} ) | \
    ( ${vaf} & ${isHet} & (countVariant() >= ${N}) ) \
    ${snpsift_optional} \
    " > ${sample_id}.germline.vcf
  >>>

  output {
    File phased_vcf = "${sample_id}.germline.vcf"
  }

  runtime {
    docker: 'docker://davelabhub/snpsift:4.3.1t--1'
    shell: ['/bin/bash', '-euo', 'pipefail']
  }
}