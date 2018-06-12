#!/bin/sh
#BSUB -J define_loci_master
#BSUB -q long
#BSUB -n 1
#BSUB -R "select[mem>2000] rusage[mem=2000] span[hosts=1]" -M2000
#BSUB -o logs/output.%J.txt
#BSUB -e logs/errorfile.%J.txt

# Source environment
source activate finemapping

# Run queue submission pipeline
snakemake \
  --rerun-incomplete \
  -s 1b.define_loci.queue_submission.Snakefile \
  --nolock \
  --local-cores 1 \
  --jobs 500 \
  --restart-times 2 \
  --cluster-config configs/cluster.json \
  --cluster "python scripts/bsub.py"
