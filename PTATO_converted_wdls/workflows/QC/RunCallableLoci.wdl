version 1.0

import "../../NextflowModules/Utils/getFilesFromDir.wdl" as Utils
import "../../NextflowModules/GATK/3.8.1/CallableLoci.wdl" as GATK
import "../../NextflowModules/Utils/AutosomalCallable.wdl" as AutosomalCallable

workflow RunCallableLoci {
  input {
    Array[File] bams
    String? callableloci_dir
    String? autosomal_callable_dir
    String out_dir
  }

  # Step 1: Extract or run CallableLoci
  if (defined(callableloci_dir)) {
    call Utils.extractCallableLociBedFromDir {
      input:
        dir = callableloci_dir
    }
    Array[Array[File]] callableloci_files = Utils.extractCallableLociBedFromDir.out
  } else {
    call GATK.CallableLoci {
      input:
        bams = bams
    }
    Array[Array[File]] callableloci_files = GATK.CallableLoci.out.map { 
      Array[File] [
        sep('/', out_dir, "QC", "CallableLoci", ~{it[0]}, basename(it[2])),
        sep('/', out_dir, "QC", "CallableLoci", ~{it[0]}, basename(it[3]))
      ]
    }
  }

  # Step 2: Extract or run AutosomalCallableLoci
  if (defined(autosomal_callable_dir)) {
    call Utils.extractAutosomalCallableLociFromDir {
      input:
        dir = autosomal_callable_dir
    }
    Array[Array[File]] autosomal_callable_files = Utils.extractAutosomalCallableLociFromDir.out
  } else {
    call AutosomalCallable.AutosomalCallableLoci {
      input:
        callableloci_files = callableloci_files
    }
    Array[Array[File]] autosomal_callable_files = AutosomalCallable.AutosomalCallableLoci.out.map { 
      Array[File] [
        sep('/', out_dir, "QC", "AutosomalCallableLoci", ~{it[0]}, basename(it[2]))
      ]
    }
  }

  # Step 3: Combine callable loci and autosomal callable files
  scatter (filePair in zip(callableloci_files, autosomal_callable_files)) {
    Array[File] combined_files = filePair.left + filePair.right
  }

  output {
    Array[Array[File]] callable_combined_files = combined_files
  }
}
