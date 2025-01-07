version 1.0

import "../../NextflowModules/Utils/getFilesFromDir.wdl" as Utils
import "../../NextflowModules/shapeit/4.2.2/shapeit.wdl" as Shapeit
import "../../NextflowModules/Utils/ABtable.wdl" as ABTable
import "../../NextflowModules/htslib/1.15/tabix.wdl" as HTSLib

workflow get_ab_tables {
  input {
    Array[File] input_germline_vcfs
    Array[File] input_somatic_vcfs
    String chroms_file_path
    String out_dir
    String? phased_vcfs_dir
  }

  # Step 1: Read chromosomes file and pair with germline VCFs
  File chroms_file = chroms_file_path
  Array[String] chroms = read_lines(chroms_file)
  Array[Array[File]] input_germline_vcfs_chroms = zip(input_germline_vcfs, chroms)

  # Step 2: Extract or create phased VCFs
  if (defined(phased_vcfs_dir)) {
    call Utils.extractPhasedVcfGzFromDir {
      input:
        dir = phased_vcfs_dir
    }
    Array[Array[File]] phased_vcfs = Utils.extractPhasedVcfGzFromDir.out
  } else {
    call Shapeit.shapeit {
      input:
        input_germline_vcfs_chroms = input_germline_vcfs_chroms
    }

    call HTSLib.tabix {
      input:
        input_vcfs = Shapeit.shapeit.out
    }

    Array[Array[File]] phased_vcfs = HTSLib.tabix.out.map { 
      Array[File] [
        sep('/', out_dir, "intermediate", "short_variants", "shapeit", ~{it[0]}, basename(it[2])),
        sep('/', out_dir, "intermediate", "short_variants", "shapeit", ~{it[0]}, basename(it[3]))
      ]
    }
  }

  # Step 3: Create AB tables
  call ABTable.createABtable {
    input:
      input_files = zip(zip(phased_vcfs, input_germline_vcfs), input_somatic_vcfs)
        .map {
          Array[File] [
            it[0][0],
            it[0][1],
            it[1]
          ]
        }
  }

  # Step 4: Merge AB tables
  call ABTable.mergeABtable {
    input:
      input_ab_tables = ABTable.createABtable.out.group_by(0, 1)
  }

  # Step 5: Copy AB tables to output directory
  Array[Array[File]] ab_tables = ABTable.mergeABtable.out.map { 
    Array[File] [
      sep('/', out_dir, "intermediate", "short_variants", "ab", ~{it[0]}, basename(it[2]))
    ]
  }

  output {
    Array[Array[File]] ab_tables
  }
}
