#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#
"""
Converts from cleaned format to GCTA format for finemapping
"""

import os
import sys
import argparse
import pandas as pd

def main():

    # Parse args
    args = parse_args()

    # Load iteratively, filter by maf
    iter_tsv = pd.read_csv(args.inf, sep="\t", header=0, iterator=True, chunksize=1000000)
    data = pd.concat([chunk[chunk['eaf'].apply(lambda x: min(x, 1-x)) >= args.maf] for chunk in iter_tsv])

    # Rename columns
    data = data.rename(columns={"rsid":"SNP",
                                "effect_allele":"A1",
                                "other_allele":"A2",
                                "eaf":"freq",
                                "beta":"b",
                                "pval":"p",
                                "n":"N"})
    data = data.loc[:, ["SNP", "A1", "A2", "freq", "b", "se", "p", "N"]]
    data.to_csv(args.outf, sep="\t", index=None, compression="gzip")

    return 0

def parse_args():
    """ Load command line args.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--inf', metavar="<file>", help=('Input file'), type=str, required=True)
    parser.add_argument('--outf', metavar="<file>", help=('Output file'), type=str, required=True)
    parser.add_argument('--maf', metavar="<float>", help=('Min MAF (default: 0.01)'), type=float, default=0.01)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
