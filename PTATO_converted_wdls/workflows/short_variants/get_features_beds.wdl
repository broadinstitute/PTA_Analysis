version 1.0

import "../../NextflowModules/bedtools/2.30.0/closest.wdl" as BedtoolsClosest
import "../../NextflowModules/bedtools/2.30.0/intersect.wdl" as BedtoolsIntersect
import "../../NextflowModules/bedtools/2.30.0/groupby.wdl" as BedtoolsGroupby

workflow closest_feature {
  input {
    Array[Array[File]] input_sample_beds
    Array[Array[File]] input_feature_beds
    String out_dir
  }

  # Combine input sample beds with feature beds
  Array[Array[Array[File]]] input_files = zip(input_sample_beds, input_feature_beds)

  # Call Bedtools Closest
  call BedtoolsClosest.closest {
    input:
      input_files = input_files
  }

  # Call Bedtools Groupby
  call BedtoolsGroupby.groupby {
    input:
      input_files = BedtoolsClosest.closest.out
  }

  # Copy features beds to output directory
  Array[Array[File]] features_beds = BedtoolsGroupby.groupby.out.map {
    Array[File] [
      sep('/', out_dir, "intermediate", "short_variants", "features", ~{it[0]}, ~{it[1]}, basename(it[2]))
    ]
  }

  output {
    Array[Array[File]] features_beds
  }
}

workflow intersect_feature {
  input {
    Array[Array[File]] input_sample_beds
    Array[Array[File]] input_feature_beds
    String out_dir
  }

  # Combine input sample beds with feature beds
  Array[Array[Array[File]]] input_files = zip(input_sample_beds, input_feature_beds)

  # Call Bedtools Intersect
  call BedtoolsIntersect.intersect {
    input:
      input_files = input_files
  }

  # Call Bedtools Groupby
  call BedtoolsGroupby.groupby {
    input:
      input_files = BedtoolsIntersect.intersect.out
  }

  # Copy features beds to output directory
  Array[Array[File]] features_beds = BedtoolsGroupby.groupby.out.map {
    Array[File] [
      sep('/', out_dir, "intermediate", "short_variants", "features", ~{it[0]}, ~{it[1]}, basename(it[2]))
    ]
  }

  output {
    Array[Array[File]] features_beds
  }
}

workflow groupby_features {
  input {
    Array[Array[File]] input_sample_beds
    Array[Array[File]] input_feature_groupby_beds
    String out_dir
  }

  # Join input sample beds with feature groupby beds
  Array[Array[Array[File]]] input_files = join_by_key(input_sample_beds, input_feature_groupby_beds, keys = [0, 1])

  # Call Bedtools IntersectAll
  call BedtoolsIntersect.intersectAll {
    input:
      input_files = input_files
  }

  # Call Bedtools GroupbyAll
  call BedtoolsGroupby.groupbyAll {
    input:
      input_files = BedtoolsIntersect.intersectAll.out
  }

  # Copy features beds to output directory
  Array[Array[File]] features_beds = BedtoolsGroupby.groupbyAll.out.map {
    Array[File] [
      sep('/', out_dir, "intermediate", "short_variants", "features", ~{it[0]}, basename(it[2]))
    ]
  }

  output {
    Array[Array[File]] features_beds
  }
}
