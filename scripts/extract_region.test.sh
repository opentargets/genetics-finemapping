#!/usr/bin/env bash
#

python extract_region.py \
  ../data/SHARE-without23andMe.gcta_format.tsv \
  ../../reference_data_download/1000Genomes_phase3/plink_format/EUR/1kg_p3.20130502.EUR.1.bim \
  rs2228145 \
  500 \
  0.01 \
  ../temp/extract_test.tsv
