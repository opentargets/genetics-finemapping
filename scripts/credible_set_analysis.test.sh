#!/usr/bin/env bash
#

# python credible_set_analysis.py \
#        --inf ../results/SHARE-without23andMe/cond_analysis/1.cond.rs1102705.conditional_sumstats.tsv \
#        --outf ../temp/1.cond.rs1102705.cred_set.tsv \
#        --prop_cases 0.6203537

python ../scripts/credible_set_analysis.py --inf ../results/20002_1226/cond_analysis/6.cond.rs117108573.conditional_sumstats.tsv --outf ../results/20002_1226/credible_set/6.cond.rs117108573.credible_sets.tsv --prop_cases 0.04857055573186538
