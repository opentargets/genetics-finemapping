#!/usr/bin/env bash
#

set -euo pipefail

# Neale
mkdir -p raw/NEALEUKB_50/UKB_50
gsutil -m cp -r "gs://genetics-portal-sumstats/gwas/genome_wide/NEALEUKB_50/UKB_50/*-NEALEUKB_50-UKB_50.tsv.gz" raw/NEALEUKB_50/UKB_50

# GTEx
mkdir -p raw/GTEX/UBERON_0000178
gsutil -m cp -r "gs://genetics-portal-raw/eqtl_gtex_v7/allpairs/Whole_Blood.allpairs.txt.gz" raw/GTEX/UBERON_0000178


echo COMPLETE
