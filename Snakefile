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

import os
import sys
import itertools
import pandas as pd
import subprocess as sp

# Args
configfile: "configs/config.yaml"
studies,  = glob_wildcards("data/{study}.txt.gz")
chroms = [str(x) for x in range(1, 23)]
# chroms = "X"
# chroms = "1"

#
# Stage 1: Define independent loci using GCTA COJO slct ------------------------
#

rule define_loci:
    """ Make targets for GCTA cojo slct. Loci need to be defined before
        finemapping stage.
    """
    input:
        expand("results/{study}/indep_loci/{chrom}.gcta_slct.jma.cojo",
               study=studies, chrom=chroms)

rule format_sumstats_for_gcta:
    """ This function is specific to my test dataset, it will need updating to
        fit the generic format
    """
    input:
        "data/{study}.txt.gz"
    output:
        "data/{study}.gcta_format.tsv"
    shell:
        'python scripts/format_for_gcta.py '
        '--inf {input} '
        '--outf {output} '
        '--snp_col RS_ID '
        '--effect_col EFFECT_ALLELE '
        '--other_col OTHER_ALLELE '
        '--freq_col 1000G_ALLELE_FREQ '
        '--beta_col BETA '
        '--se_col SE '
        '--p_col PVALUE '
        '--n_col N '
        '--sep " "'

rule cojo_indep_loci:
    """ Use GCTA COJO slct to identify independent loci
    """
    input:
        "data/{study}.gcta_format.tsv"
    output:
        "results/{study}/indep_loci/{chrom}.gcta_slct.jma.cojo"
    params:
        outpref="results/{study}/indep_loci/{chrom}.gcta_slct",
        ref=config["ref_panel"],
        colin=config["cojo_colin"],
        p=config["cojo_p"],
        wind=config["cojo_wind"],
        maf=config["cojo_maf"]
    shell:
        'gcta64 --bfile {params.ref} '
        '--chr {wildcards.chrom} '
        '--maf {params.maf} '
        '--cojo-p {params.p} '
        '--cojo-wind {params.wind} '
        '--cojo-collinear {params.colin} '
        '--cojo-file {input} '
        '--cojo-slct '
        '--out {params.outpref}'

#
# Stage 2: Credible set analysis for loci defined in stage 1 -------------------
#

def make_study_loci_dict(in_format="results/{study}/indep_loci/{chrom}.gcta_slct.jma.cojo"):
    """ Reads the gcta jma (indepent loci file) to make a dictionary of loci
        for each study.
    """
    d = {}
    for study, chrom in itertools.product(studies, chroms):
        if study not in d: d[study] = {}
        if chrom not in d[study]: d[study][chrom] = []
        in_file = in_format.format(study=study, chrom=chrom)
        if os.path.exists(in_file):
            with open(in_file, "r") as in_jma:
                in_jma.readline() # Skip header
                for line in in_jma:
                    index_snp = line.rstrip().split("\t")[1]
                    d[study][chrom].append(index_snp)
    return d

# Get dictionary of independent loci per study per chrom
study_loci = make_study_loci_dict(in_format="results/{study}/indep_loci/{chrom}.gcta_slct.jma.cojo")

def make_finemap_targets(study_loci, out_format):
    """ Reads GCTA cojo jma files to get index_snps for each study /chromosome.
        Uses these to produce target files for finemapping
    """
    outfiles = []
    for study in study_loci:
        for chrom in study_loci[study]:
            for index_snp in study_loci[study][chrom]:
                outfiles.append(out_format.format(study=study,
                                                  chrom=chrom,
                                                  index_snp=index_snp))
    return outfiles

rule finemap_loci:
    """ Make targets for finemapping.
    """
    input:
        # make_finemap_targets(
        #     study_loci,
        #     out_format="results/{study}/credible_set/{chrom}.cond.{index_snp}.credible_sets.tsv")
        expand("results/{study}/stats/histogram.credible_set.stats.txt", study=studies)

rule identify_multi_sig_loci:
    """ Clusters the independent loci to identify SNPs within N KB of each other
    """
    input:
        "results/{study}/indep_loci/{chrom}.gcta_slct.jma.cojo"
    output:
        "results/{study}/cond_analysis/{chrom}.cond.config.txt"
    params:
        wind=config["cond_multisig_wind"]
    shell:
        'python scripts/find_multiple_signals.py {input} {params.wind} {output}'

rule write_conditional_snplist:
    """ For each index snp, write a file containing other SNPs on which the
        summary stats should be conditioned using GCTA COJO
    """
    input:
        "results/{study}/cond_analysis/{chrom}.cond.config.txt"
    output:
        "results/{study}/cond_analysis/{chrom}.cond.{index_snp}.snplist.txt"
    shell:
        'python scripts/make_conditional_snp_file.py {input} {output} {wildcards.index_snp}'

def ismultisig(inf):
    """ Determine if locus needs conditional analysis to be run (cond list will
        have >0 lines)
    """
    c = 0
    with open(inf, "r") as in_h:
        for line in in_h:
            c += 1
    return c > 0

rule run_conditional_analysis:
    """ Run conditional analysis using GCTA COJO. Conditional analysis on needs
        to be run if there are any snps contained in the *.snplist.txt file.
        Else, extract uncorrected summary stats.
    """
    input:
        sumstats="data/{study}.gcta_format.tsv",
        condlist="results/{study}/cond_analysis/{chrom}.cond.{index_snp}.snplist.txt"
    output:
        "results/{study}/cond_analysis/{chrom}.cond.{index_snp}.conditional_sumstats.tsv"
    params:
        ref=config["ref_panel"],
        maf=config["cond_maf"],
        window=config["cond_extract_wind"],
        prefix="results/{study}/cond_analysis/{chrom}.cond.{index_snp}"
    run:

        # If the locus has multiple independent signals, run condition analysis
        if ismultisig(input["condlist"]):
            cmd =  ["gcta64 --bfile {0}".format(params["ref"]),
                    "--chr {0}".format(wildcards["chrom"]),
                    "--maf {0}".format(params["maf"]),
                    "--extract-region-snp {0} {1}".format(wildcards["index_snp"], params["window"]),
                    "--cojo-file {0}".format(input["sumstats"]),
                    "--cojo-cond {0}".format(input["condlist"]),
                    "--out {0}".format(params["prefix"])]
            sp.run(" ".join(cmd), shell=True)
            # Rename from *.cma.cojo to *.conditional_sumstats.tsv
            os.rename("{0}.cma.cojo".format(params["prefix"]), str(output))

        # Otherwise, output sumstats in the same format as conditional analysis
        else:
            cmd = ["python scripts/extract_region.py",
                   input["sumstats"],
                   "{0}.bim".format(params["ref"]),
                   wildcards["index_snp"],
                   params["window"],
                   params["maf"],
                   output]
            sp.run(" ".join(str(x) for x in cmd), shell=True)

rule credible_set_analysis:
    """ Runs credible set analysis for each indepenedent locus using summary
        stats that are conditional on other nearby index snps
    """
    input:
        "results/{study}/cond_analysis/{chrom}.cond.{index_snp}.conditional_sumstats.tsv"
    output:
        "results/{study}/credible_set/{chrom}.cond.{index_snp}.credible_sets.tsv"
    params:
        prop_cases=0.6203537 # TODO get actual prop_cases from somewhere
    shell:
        "python scripts/credible_set_analysis.py "
        "--inf {input} "
        "--outf {output} "
        "--prop_cases {params.prop_cases}"

def list_merge_input_files(study, informat="results/{study}/credible_set/{chrom}.cond.{index_snp}.credible_sets.tsv"):
    for chrom in study_loci[study]:
        for index_snp in study_loci[study][chrom]:
            yield informat.format(study=study, chrom=chrom, index_snp=index_snp)

rule merge_to_bed_per_study:
    """ Merges results of credible set analysis into a single bed file per
        study
    """
    input:
        lambda wildcards: list(list_merge_input_files(wildcards.study, informat="results/{study}/credible_set/{chrom}.cond.{index_snp}.credible_sets.tsv"))
    output:
        "results/{study}/credible_set/merged.{credset}.credible_set.bed"
    shell:
        "python scripts/credset_to_bed.py "
        "--infiles {input} "
        "--outbed {output} "
        "--credset {wildcards.credset}"

rule credset_stats:
    """ Calc average set size, etc.. and make plots
    """
    input:
        "results/{study}/credible_set/merged.95.credible_set.bed",
        "results/{study}/credible_set/merged.99.credible_set.bed"
    output:
        # plot="results/{study}/stats/histogram.credible_set.png",
        stat="results/{study}/stats/histogram.credible_set.stats.txt"
    shell:
        "python scripts/plot_histogram.py --infiles {input} "
        "--out {output.stat} "
        # "--plot {output.plot} "
        "--names 95% 99%"
