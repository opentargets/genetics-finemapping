Finemapping
===========

Finemapping pipeline for summary statistics. Currently only implements credible set analysis.

TODO:
  - Write documentation
  - Get proportion of cases for each study
  - When using better reference, RAM requirements will need increasing
  - When using non-UKB sum stats, proportion of cases needs to be calculated from somewhere

Issues:
  - Chromosome X currently isn't working (error code 134), there are no X chrom SNPs in input
  - Stage 1 currently fails if no SNPs are selected

Questions:
  - What do we want the threshold for defining independent loci to be? Currently set at 0.8.

Requirements
  - GCTA v1.91.1b (there is a bug in v1.91.2b which causes GCTA slct to produce unexpected results)

### Usage

```
# Install dependencies into isolated environment
conda env create -n finemapping --file environment.yaml

# Activate environment
source activate finemapping

# Alter configuration file
nano configs/config.yaml

# Execute workflow (locally)
snakemake -p define_loci
snakemake -p finemap_loci

# Execute workflow (on Sanger farm)
bsub < bsub_wrapper.sh
```
