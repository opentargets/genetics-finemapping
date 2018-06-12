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
studies = [config["study"].replace("f.", "")]

# Define independent loci using GCTA COJO slct ---------------------------------

rule all:
    """ Make targets for GCTA cojo slct. Loci need to be defined
        before finemapping stage.

    """
    input:
        expand("results/{study}/indep_loci_cojo/{chrom}.gcta_slct.jma.cojo",
               study=studies, chrom=config["chroms"])

rule nealeUKB_to_gcta:
    """ Converts from Neale UKB cleaned format to GCTA format for finemapping.
        Takes ~5 mins.
    """
    input:
        config["in_sumstats"]
    output:
        "temp/{study}.gcta_format.tsv.gz"
    params:
        maf=config["cojo_maf"]
    shell:
        "python scripts/nealeUKB_to_gcta.py --inf {input} "
        "--outf {output} "
        "--maf {params.maf}"

rule unzip_sumstats:
    """ Locale rule. Helper rule to unzip the input for GCTA cojo """
    input:
        "temp/{study}.gcta_format.tsv.gz"
    output:
        temp("temp/{study}.gcta_format.tsv")
    shell:
        "zcat < {input} > {output}"

rule cojo_indep_loci:
    """ Use GCTA COJO slct to identify independent loci.
    Issues:
      - If no SNPs are selected, it exits with error code 0 but doesn't create
        an output file. This is problematic for snakemake.
      - If SNPs for a given chrom (e.g. X) don't exist in the input, GCTA
        returns a 134 error (program aborted) and doesn't create an output
    """
    input:
        "temp/{study}.gcta_format.tsv"
    output:
        "results/{study}/indep_loci_cojo/{chrom}.gcta_slct.jma.cojo"
    params:
        outpref="results/{study}/indep_loci_cojo/{chrom}.gcta_slct",
        ref=config["ref_panel"],
        colin=config["cojo_colin"],
        p=config["cojo_p"],
        wind=config["cojo_wind"],
        maf=config["cojo_maf"],
        exl_mhc=config["exclude_MHC"]
    run:
        # Run GCTA slct
        # cmd = ['gcta64 --bfile {0}'.format(params["ref"]),
        cmd = ['gcta64_1.91.4b --bfile {0}'.format(params["ref"]),
               '--chr {0}'.format(wildcards["chrom"]),
               '--maf {0}'.format(params["maf"]),
               '--cojo-p {0}'.format(params["p"]),
               '--cojo-wind {0}'.format(params["wind"]),
               '--cojo-collinear {0}'.format(params["colin"]),
               '--cojo-file {0}'.format(input[0]),
               '--cojo-slct',
               '--out {0}'.format(params["outpref"])]
        # Exclude MHC
        if params["exl_mhc"] and wildcards["chrom"] == "6":
            cmd.append("--exclude-region-bp 6 30963074 2486")
        print(" ".join(cmd))
        cp = sp.run(" ".join(cmd), shell=True)

        # If no independent loci found, create empty output file
        header = ["Chr", "SNP", "bp", "refA", "freq", "b", "se", "p", "n", "freq_geno", "bJ", "bJ_se", "pJ", "LD_r"]
        if cp.returncode == 0 and not os.path.exists(output[0]):
            with open(output[0], "w") as out_h:
                out_h.write("\t".join(header) + "\n")
            print("\nWarning: No independent loci found when creating {0}\n".format(output[0]))
        elif cp.returncode != 0:
            sys.exit("Error: GCTA finished with error code {0} when creating {1}".format(cp.returncode, output[0]))
