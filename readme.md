Finemapping pipleine
====================

Finemapping pipeline for Open Targets Genetics. The pipeline is currently bespoke to the Neale lab UK Biobank summary statistics (version 1) and must be run on a cluster using `bsub`.

A generalised finemapping/colocalisation pipeline is under development ([analysis plan](https://docs.google.com/document/d/1m2XFvovzXtFKoH9J-aZriH69k54Wa5cYwQJEh7VrJAs/edit?usp=sharing)).

### Contents
- [Requirements](#requirements)
- [Usage](#usage)
- [Methods](#methods)

### Requirements
- GCTA >= v1.91.3 (there is a bug in v1.91.2b which causes GCTA slct to produce unexpected results)
- `conda`

### Usage

```
# Install dependencies into isolated environment
conda env create -n finemapping --file environment.yaml

# Activate environment
source activate finemapping

# Alter configuration file
nano configs/config.yaml

# Execute steps sequentially (on cluster using BSUB)
bsub < 1a.define_loci.bsub_wrapper.sh
bsub < 2a.finemap_loci_multisig.bsub_wrapper.sh
bsub -J collate_results -n 4 -q normal -R "select[mem>16000] rusage[mem=16000] span[hosts=1]" -M16000 snakemake -s 3.collate_finemapping_results.Snakefile --cores 4

```

### Methods

#### Overview

![workflow-overview](overview_figure.png)
