#!/usr/bin/env snakemake
#
# Ed Mountjoy
#

import os
import sys

configfile: "configs/config.yaml"

# Get list of studies from GCTA cojo files
studies, _ = glob_wildcards("results/{study}/indep_loci_cojo/{chrom}.gcta_slct.jma.cojo")
studies = list(set(studies))
#studies = ["50"]

# Define rules that should be run locally by the master script, rather than distributed
localrules: all

rule all:
    """ Make targets for finemapping.
    """
    input:
        expand("results/{study}/credible_set/merged.95.credible_set.bed", study=studies)

# Create target files using snakemake pipeline (1 per study)
rule run_finemapping_pipeline:
    """ Runs the Snakemake finemapping pipeline
    """
    output:
        "results/{study}/credible_set/merged.95.credible_set.bed"
    params:
        # There is a bug in snakemake that means underscores are not
        # preserved when using --config:
        # https://bitbucket.org/snakemake/snakemake/issues/795/underscore-and-number-from-string
        # Therefore replacing prefixing with "f."
        study = lambda wildcards: "f." + wildcards["study"]
    threads: config["n_cores_finemap"]
    resources:
        mem=4000 * config["n_cores_finemap"],
        runtime=lambda wildcards, attempt: 60 * 12 if attempt == 1 else 60 * 48
    shell:
        "snakemake -s 2c.finemap_loci_multisig.pipeline.Snakefile "
        "--nolock "
        "--cores {threads} "
        "--config study={params.study}"
