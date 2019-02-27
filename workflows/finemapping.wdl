# Workflow combining tasks together
workflow finemapping {
  File inputManifestFile
  Array[Array[String]] inputManifest = read_tsv(inputManifestFile)
  scatter (study in inputManifest) {
    call finemap_single_study {
      input:
        type=study[0],
        in_pq=study[1],
        in_ld=study[2],
        study_id=study[3],
        phenotype_id=study[4],
        biofeature=study[5],
        chrom=study[6],
        out_top_loci=study[7],
        out_credset=study[8],
        out_log=study[9],
        out_tmpdir=study[10],
        method=study[11],
        pval_threshold=study[12]
    }
  }
  call concat_json as concat_toploci {
    input:
      in=finemap_single_study.toploci_res
  }
  call concat_json as concat_credset {
    input:
      in=finemap_single_study.credset_res
  }
}

# Task to run the finemapping script
task finemap_single_study {
  
  # Arguments specified in workflow.config.json
  String script
  File config_file

  # Arguments specified in input_files.config.tsv
  String type
  String in_pq
  String in_ld
  String study_id
  String phenotype_id
  String biofeature
  String chrom
  String out_top_loci
  String out_credset
  String out_log
  String out_tmpdir
  String method
  Float pval_threshold

  command {
    python ${script} \
      --pq ${in_pq} \
      --ld ${in_ld} \
      --config_file ${config_file} \
      --type ${type} \
      --study_id ${study_id} \
      --phenotype_id ${phenotype_id} \
      --biofeature ${biofeature} \
      --chrom ${chrom} \
      --toploci ${out_top_loci} \
      --credset ${out_credset} \
      --log ${out_log} \
      --tmpdir ${out_tmpdir} \
      --method ${method} \
      --pval_threshold ${pval_threshold}
  }
  output {
    File toploci_res = "${out_top_loci}"
    File credset_res = "${out_credset}"
  }
}

# Concatenate jsons together
task concat_json {
  Array[File] in
  String out
  command {
    zcat < ${sep=" " in} | gzip -c > ${out}
  }
  output {
    File result = "${out}"
  }
}

# # Task to load parquets into memory and write to a single file
# task concat_parquets {
#   String script
#   Array[File] in
#   String out
#   command {
#     python ${script} \
#       --in_parquets ${sep=" " in} \
#       --out ${out} \
#   }
#   output {
#     File result = "${out}"
#   }
# }
