#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#
"""
Format summary stats for use with GCTA's cojo. Requires:
  SNP A1 A2 freq b se p N
  where:
    A1 is the effect allele
    freq is wrt A1
    b and se are log(OR) for case-control studies
"""

import sys
import gzip
import argparse
import pandas as pd

def main():

    # Load args
    args = parse_args()

    # Load data to pandas
    data = pd.read_csv(args.inf, header=0, sep=args.sep)

    # Rename columns
    data = data.rename(columns={"SNP":"CHR_POS"})
    data = data.rename(columns={args.snp_col:"SNP",
                                args.effect_col:"A1",
                                args.other_col:"A2",
                                args.freq_col:"freq",
                                args.beta_col:"b",
                                args.se_col:"se",
                                args.p_col:"p",
                                args.n_col:"N"})

    # Make alleles upper case
    data.loc[:, "A1"] = data.loc[:, "A1"].apply(lambda x: x.upper())
    data.loc[:, "A2"] = data.loc[:, "A2"].apply(lambda x: x.upper())

    # Flip allele frequencies to match A1
    to_fix = data.A1 != data.loc[:, "1000G_ALLELE"]
    data.loc[to_fix, "freq"] = data.loc[to_fix, "freq"].apply(lambda x: 1 - x)

    data.loc[:, ["SNP", "A1", "A2", "freq", "b", "se", "p", "N"]].to_csv(args.outf, sep="\t", index=False)

    return 0

def parse_args():
    """ Load command line args.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--inf', metavar="<file>", help=('GWAS summary statistics file'), type=str, required=True)
    parser.add_argument('--outf', metavar="<file>", help=("Formatted summary statistics file"), type=str, required=True)
    parser.add_argument('--snp_col', metavar="<file>", help=('SNP col'), type=str, required=True)
    parser.add_argument('--effect_col', metavar="<file>", help=('Effect allele col'), type=str, required=True)
    parser.add_argument('--other_col', metavar="<file>", help=('Effect allele col'), type=str, required=True)
    parser.add_argument('--freq_col', metavar="<file>", help=('Freq of effect allele col'), type=str, required=True)
    parser.add_argument('--beta_col', metavar="<file>", help=('Beta col'), type=str, required=True)
    parser.add_argument('--se_col', metavar="<file>", help=('SE col'), type=str, required=True)
    parser.add_argument('--p_col', metavar="<file>", help=('P col'), type=str, required=True)
    parser.add_argument('--n_col', metavar="<file>", help=('N col'), type=str, required=False)
    parser.add_argument('--n', metavar="<file>", help=('Sample N over-ride'), type=int, required=False)
    parser.add_argument('--sep', metavar="<file>", help=('Column sep (default: tab)'), type=str, default="\t")

    args = parser.parse_args()
    if not args.n_col and not args.n:
        sys.exit("Error: either --n_col or --n are required")

    return args


if __name__ == '__main__':

    main()
