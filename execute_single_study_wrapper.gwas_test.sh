#!/usr/bin/env bash
#

set -euo pipefail

# Run example gwas study
python finemapping/single_study.wrapper.py \
  --pq 'example_data/sumstats/gwas/GCST004132_cr.parquet' \
  --ld '/Users/em21/Projects/reference_data/uk10k_2019Feb/3_liftover_to_GRCh38/output/{chrom}.ALSPAC_TWINSUK.maf01.beagle.csq.shapeit.20131101' \
  --config_file 'configs/analysis.config.yaml' \
  --study_id 'GCST004132_cr' \
  --phenotype_id 'None' \
  --bio_feature 'None' \
  --type 'gwas' \
  --chrom '22' \
  --method 'conditional' \
  --pval_threshold 5e-8 \
  --toploci 'output/study_id=GCST004132_cr/phenotype_id=/bio_feature=/chrom=22/top_loci.json.gz' \
  --credset 'output/study_id=GCST004132_cr/phenotype_id=/bio_feature=/chrom=22/credible_set.json.gz' \
  --tmpdir 'tmp/study_id=GCST004132_cr/phenotype_id=/bio_feature=/chrom=22/' \
  --log 'logs/study_id=GCST004132_cr/phenotype_id=/bio_feature=/chrom=22/logfile.txt'
