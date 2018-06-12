#!/usr/bin/env bash
#

python format_for_gcta.py \
  --inf ../data/sumstats/SHARE-without23andMe.LDSCORE-GC.SE-META.v0.gz \
  --outf ../temp/SHARE-without23andMe.gcta_format.tsv \
  --snp_col RS_ID \
  --effect_col EFFECT_ALLELE \
  --other_col OTHER_ALLELE \
  --freq_col 1000G_ALLELE_FREQ \
  --beta_col BETA \
  --se_col SE \
  --p_col PVALUE \
  --n_col N \
  --sep " "
