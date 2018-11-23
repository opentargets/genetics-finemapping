#!/usr/bin/env bash
#

set -euo pipefail

# Run example gwas study
python finemapping/single_study.wrapper.py \
  --pq 'input/molecular_qtl/GTEX7' \
  --ld 'input/ld/EUR.{chrom}.1000Gp3.20130502' \
  --config_file 'configs/analysis.config.yaml' \
  --study_id 'GTEX7' \
  --cell_id 'UBERON_0000178' \
  --group_id 'ENSG00000258289' \
  --trait_id 'eqtl' \
  --chrom '14' \
  --method 'conditional' \
  --toploci 'output/study_id=GTEX7/cell_id=UBERON_0000178/group_id=ENSG00000258289/trait_id=eqtl/chrom=14/top_loci.parquet' \
  --credset 'output/study_id=GTEX7/cell_id=UBERON_0000178/group_id=ENSG00000258289/trait_id=eqtl/chrom=14/credible_set.parquet' \
  --tmpdir 'tmp/study_id=GTEX7/cell_id=UBERON_0000178/group_id=ENSG00000258289/trait_id=eqtl/chrom=14/' \
  --log 'logs/study_id=GTEX7/cell_id=UBERON_0000178/group_id=ENSG00000258289/trait_id=eqtl/chrom=14/logfile.txt'

echo COMPLETE
