Finemapping
===========

Implements a generalised credible set analysis

TODO:
  - Allow calculating of approximate Bayes factors for quantitative traits, currently limited to case/control studies
  - Write documentation

# install dependencies into isolated environment
conda env create -n finemapping --file environment.yaml

# activate environment
source activate finemapping

# execute workflow
bash run_finemap_pipeline.sh
