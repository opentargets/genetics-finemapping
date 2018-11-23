#!/usr/bin/env bash
#

set -euo pipefail

# Run example gwas study
python finemapping/single_study.wrapper.py \
  --pq 'input/gwas/NEALEUKB_50' \
  --ld 'input/ld/EUR.{chrom}.1000Gp3.20130502' \
  --config_file 'configs/analysis.config.yaml' \
  --study_id 'NEALEUKB_50' \
  --trait_id 'UKB_50' \
  --chrom '22' \
  --method 'conditional' \
  --toploci 'output/study_id=NEALEUKB_50/cell_id=/group_id=/trait_id=UKB_50/chrom=22/top_loci.parquet' \
  --credset 'output/study_id=NEALEUKB_50/cell_id=/group_id=/trait_id=UKB_50/chrom=22/credible_set.parquet' \
  --tmpdir 'tmp/study_id=NEALEUKB_50/cell_id=/group_id=/trait_id=UKB_50/chrom=22/' \
  --log 'logs/study_id=NEALEUKB_50/cell_id=/group_id=/trait_id=UKB_50/chrom=22/logfile.txt'

  # Not required
  # --cell_id '' \
  # --group_id '' \
  # --trait_id '' \

echo COMPLETE
