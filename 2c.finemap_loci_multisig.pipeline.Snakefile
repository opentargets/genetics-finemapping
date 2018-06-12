#!/usr/bin/env snakemake
#
# Ed Mountjoy
#
# Snakemake pipeline for credible set analysis using summary statistics.
#

import os
import sys
import itertools
import pandas as pd
import subprocess as sp

# Args
configfile: "configs/config.yaml"
studies = [config["study"].replace("f.", "")]

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

def make_study_loci_dict_cojo(in_format="results/{study}/indep_loci_cojo/{chrom}.gcta_slct.jma.cojo"):
    """ Reads the gcta jma (indepent loci file) to make a dictionary of loci
        for each study.
    Returns:
        dict {study:{chrom: [index_snp1, ...]}}
    """
    d = {}
    for study, chrom in itertools.product(studies, config["chroms"]):
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
study_loci_cojo = make_study_loci_dict_cojo(in_format="results/{study}/indep_loci_cojo/{chrom}.gcta_slct.jma.cojo")

rule finemap_loci_multisig:
    """ Make targets for finemapping.
    """
    input:
        expand("results/{study}/credible_set/merged.{credset}.credible_set.bed", study=studies, credset=[95, 99])

rule nealeUKB_to_gcta:
    """ Converts from Neale UKB cleaned format to GCTA format for finemapping.
        Takes ~5 mins.
    """
    input:
        config["in_sumstats"]
    output:
        temp("temp/{study}.gcta_format.tsv.gz")
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

rule identify_multi_sig_loci:
    """ Local rule. Clusters the independent loci to identify SNPs within N KB
        of each other
    """
    input:
        "results/{study}/indep_loci_cojo/{chrom}.gcta_slct.jma.cojo"
    output:
        "results/{study}/cond_analysis/{chrom}.cond.config.txt"
    params:
        wind=config["cond_multisig_wind"]
    shell:
        'python scripts/find_multiple_signals.py {input} {params.wind} {output}'

rule write_conditional_snplist:
    """ Local rule. For each index snp, write a file containing other SNPs on
        which the summary stats should be conditioned using GCTA COJO
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
        sumstats="temp/{study}.gcta_format.tsv",
        condlist="results/{study}/cond_analysis/{chrom}.cond.{index_snp}.snplist.txt"
    output:
        temp("results/{study}/cond_analysis/{chrom}.cond.{index_snp}.conditional_sumstats.tsv")
    params:
        ref=config["ref_panel"],
        maf=config["cond_maf"],
        window=config["cond_extract_wind"],
        prefix="results/{study}/cond_analysis/{chrom}.cond.{index_snp}"
    resources:
        mem=4000, # Guess
        runtime=lambda wildcards, attempt: attempt * 10
    run:

        # If the locus has multiple independent signals, run condition analysis
        if ismultisig(input["condlist"]):
            # cmd =  ["gcta64 --bfile {0}".format(params["ref"]),
            cmd =  ["gcta64_1.91.4b --bfile {0}".format(params["ref"]),
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
    """ Local rule. Runs credible set analysis for each indepenedent locus using
        summary stats that are conditional on other nearby index snps
    """
    input:
        "results/{study}/cond_analysis/{chrom}.cond.{index_snp}.conditional_sumstats.tsv"
    output:
        "results/{study}/credible_set/{chrom}.cond.{index_snp}.credible_sets.tsv"
    params:
        # prop_cases=0.6203537 # DEBUG
        prop_cases=lambda wildcards: pheno_config[wildcards.study]["case_prop"]
    run:
        # Build command
        cmd = ["python scripts/credible_set_analysis.py",
               "--inf {0}".format(input[0]),
               "--outf {0}".format(output[0])]
        # Add proportion of cases only if CC
        if pd.notnull(params["prop_cases"]):
            cmd.append("--prop_cases {0}".format(params["prop_cases"]))
        # Execute
        print(" ".join(cmd)) # DEBUG
        sp.run(" ".join(cmd), shell=True)

def list_merge_input_files(study, informat="results/{study}/credible_set/{chrom}.cond.{index_snp}.credible_sets.tsv"):
    for chrom in study_loci_cojo[study]:
        for index_snp in study_loci_cojo[study][chrom]:
            yield informat.format(study=study, chrom=chrom, index_snp=index_snp)

rule merge_to_bed_per_study:
    """ Merges results of credible set analysis into a single bed file per
        study. If there are no credible sets, will create a file with just a
        header.
    """
    input:
        lambda wildcards: list(list_merge_input_files(wildcards.study, informat="results/{study}/credible_set/{chrom}.cond.{index_snp}.credible_sets.tsv"))
    output:
        "results/{study}/credible_set/merged.{credset}.credible_set.bed"
    run:
        if len(input) > 0:
            cmd = ["python scripts/credset_to_bed.py ",
                   "--infiles {0} ".format(" ".join(input)),
                   "--outbed {0} ".format(output[0]),
                   "--credset {0}".format(wildcards["credset"])]
            print(" ".join(cmd))
            sp.run(" ".join(cmd), shell=True)
        else:
            with open(output[0], "w") as out_h:
                # Write header only
                header = ["chrom", "start", "end", "index", "size", "set"]
                out_h.write("\t".join(header) + "\n")



# rule merge_to_bed_per_study:
#     """ Local rule. Merges results of credible set analysis into a single bed
#         file per study
#     """
#     input:
#         lambda wildcards: list(list_merge_input_files(wildcards.study, informat="results/{study}/credible_set/{chrom}.cond.{index_snp}.credible_sets.tsv"))
#     output:
#         "results/{study}/credible_set/merged.{credset}.credible_set.bed"
#     shell:
#         "python scripts/credset_to_bed.py "
#         "--infiles {input} "
#         "--outbed {output} "
#         "--credset {wildcards.credset}"
