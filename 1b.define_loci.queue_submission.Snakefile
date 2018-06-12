#!/usr/bin/env snakemake
#
# Ed Mountjoy
#
# Snakemake pipeline for credible set analysis using summary statistics.
#
# Usage:
#   snakemake define_loci
#   snakemake finemap_loci
#
# Script is currently configured to process Neale UKB sumstats

import os
import sys
import itertools
import pandas as pd
import subprocess as sp

# Args
configfile: "configs/config.yaml"
studies,  = glob_wildcards(config["in_sumstats"])
# studies = ["20002_1226"] # Hyothyroidism

# Define rules that should be run locally by the master script, rather than distributed
localrules: all

# Get proportion of cases from the Neale pheno summary file
def parse_pheno_config(inf):
    """ Parses the Neale pheno summary (phenosummary_final_11898_18597.tsv)
        Need to calculate:
          - Proportion of cases (for finemap ABF calculation)
    """
    # Load
    pheno_df = pd.read_csv(inf, sep="\t", header=0)

    # Calc case prop
    pheno_df["case_prop"] = (pheno_df["N.cases"] / (pheno_df["N.cases"] + pheno_df["N.controls"]))

    # Make dict
    dict_df = pheno_df.loc[:, ["case_prop", "N.cases"]]
    dict_df.index = pheno_df["Field.code"]
    return dict_df.to_dict("index")
pheno_config = parse_pheno_config(config["pheno_file"])

# Get list of studies where N.cases > min_cases
studies_todo = [study for study in studies if
                  pheno_config[study]["N.cases"] >= config["min_cases"] or
                  pd.isnull(pheno_config[study]["N.cases"])]
studies = studies_todo

# Define targets
rule all:
    """ Local rule. Make targets for GCTA cojo slct. Loci need to be defined
        before finemapping stage.
    """
    input:
        expand("results/{study}/indep_loci_cojo/{chrom}.gcta_slct.jma.cojo",
               study=studies, chrom=config["chroms"])

# Create target files using snakemake pipeline (1 per study)
rule run_define_loci_pipeline:
    """ Runs the Snakemake pipeline
    """
    output:
        expand("results/{{study}}/indep_loci_cojo/{chrom}.gcta_slct.jma.cojo",
               chrom=config["chroms"])
    params:
        # There is a bug in snakemake that means underscores are not
        # preserved when using --config:
        # https://bitbucket.org/snakemake/snakemake/issues/795/underscore-and-number-from-string
        # Therefore replacing prefixing with "f."
        study = lambda wildcards: "f." + wildcards["study"]
    threads: config["n_cores_cojo"]
    resources:
        mem=4000 * config["n_cores_cojo"],
        runtime=lambda wildcards, attempt: 60 * 12 if attempt == 1 else 60 * 48
    shell:
        "snakemake -s 1c.define_loci.pipeline.Snakefile "
        "--nolock "
        "--cores {threads} "
        "--config study={params.study}"
