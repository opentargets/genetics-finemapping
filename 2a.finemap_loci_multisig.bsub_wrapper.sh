#!/bin/sh
#BSUB -J finemap_master
#BSUB -q long
#BSUB -n 1
#BSUB -R "select[mem>2000] rusage[mem=2000] span[hosts=1]" -M2000
#BSUB -o logs/output.%J.txt
#BSUB -e logs/errorfile.%J.txt

# Source environment
source activate finemapping

# Run queue submission pipeline
snakemake \
  -s 2b.finemap_loci_multisig.queue_submission.Snakefile \
  --nolock \
  --local-cores 1 \
  --jobs 500 \
  --restart-times 2 \
  --cluster-config configs/cluster.json \
  --cluster "python scripts/bsub.py"
