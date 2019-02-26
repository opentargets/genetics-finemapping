#!/usr/bin/env bash
#

set -euo pipefail

# Run example molecular trait study
python finemapping/single_study.wrapper.py \
  --pq 'example_data/sumstats/molecular_trait/CEDAR.parquet' \
  --ld '/Users/em21/Projects/reference_data/uk10k_2019Feb/3_liftover_to_GRCh38/output/{chrom}.ALSPAC_TWINSUK.maf01.beagle.csq.shapeit.20131101' \
  --config_file 'configs/analysis.config.yaml' \
  --study_id 'CEDAR' \
  --phenotype_id 'ILMN_1690982' \
  --bio_feature 'MONOCYTE_CD14' \
  --type 'eqtl' \
  --chrom '22' \
  --method 'conditional' \
  --toploci 'output/study_id=CEDAR/phenotype_id=ILMN_1690982/bio_feature=MONOCYTE_CD14/chrom=22/top_loci.json.gz' \
  --credset 'output/study_id=CEDAR/phenotype_id=ILMN_1690982/bio_feature=MONOCYTE_CD14/chrom=22/credible_set.json.gz' \
  --tmpdir 'tmp/study_id=CEDAR/phenotype_id=ILMN_1690982/bio_feature=MONOCYTE_CD14/chrom=22/' \
  --log 'logs/study_id=CEDAR/phenotype_id=ILMN_1690982/bio_feature=MONOCYTE_CD14/chrom=22/logfile.txt'

