Finemapping
===========

Finemapping pipeline for summary statistics. Currently only implements credible set analysis.

TODO:
  - Write documentation
  - Get proportion of cases for each study
  - Combine credible set results into a single file per study

### Usage

```
# Install dependencies into isolated environment
conda env create -n finemapping --file environment.yaml

# Activate environment
source activate finemapping

# Alter configuration file
nano configs/config.yaml

# Execute workflow
bash run_finemap_pipeline.sh
```
