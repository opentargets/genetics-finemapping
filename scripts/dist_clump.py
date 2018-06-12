#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#
"""
Performs distance based clumping. Input files must have these cols:
    snpid, chrom, pos, pval, eaf

If using --exclMHC, chrom and pos should be on assembly GRCh37. Will also
exclude long inversion on chrom 17:
    - 6:28,477,796-33,448,353
    - 17:43,384,863-44,913,631

"""

import os
import sys
import gzip
import argparse
import pandas as pd

def main():

    # Parse args
    args = parse_args()

    # Load iteratively
    iter_tsv = pd.read_csv(args.inf, sep="\t", header=0, iterator=True,
                           chunksize=1000000)
    data = pd.concat([chunk[chunk['pval'] < args.pval] for chunk in iter_tsv])
    data["chrom"] = data["chrom"].astype(str)
    if args.chrom:
        data = data.loc[data["chrom"] == args.chrom, :]
    data["maf"] = data["eaf"].apply(lambda x: min(x, 1-x))

    # Exclude MHC
    if args.exclMHC:
        to_exclude_MHC = ( (data["chrom"] == "6") &
                           (data["pos"] >= 28477796) &
                           (data["pos"] <= 33448353) )
        to_exclude_17  = ( (data["chrom"] == "17") &
                           (data["pos"] >= 43384863) &
                           (data["pos"] <= 44913631) )
        to_exclude = (to_exclude_MHC | to_exclude_17)
        data = data.loc[~to_exclude, :]

    # Process combined variants
    if args.out_combined:
        clump_combined = distance_clumping(data.copy(), args.window)
        clump_combined_index = clump_combined.loc[~clump_combined["cluster"].duplicated(keep="first"), :]
        clump_combined_index.to_csv(args.out_combined, sep="\t", index=None)

    # Process common variants
    if args.out_common:
        clump_common = distance_clumping(data.loc[data["maf"] >= 0.01, :].copy(), args.window)
        clump_common_index = clump_common.loc[~clump_common["cluster"].duplicated(keep="first"), :]
        clump_common_index.to_csv(args.out_common, sep="\t", index=None)

    # Process rare variants
    if args.out_rare:
        clump_rare = distance_clumping(data.loc[data["maf"] < 0.01, :].copy(), args.window)
        clump_rare_index = clump_rare.loc[~clump_rare["cluster"].duplicated(keep="first"), :]
        clump_rare_index.to_csv(args.out_rare, sep="\t", index=None)

    return 0

def distance_clumping(df, dist=500):
    """ Does distance based clumping.
    Args:
        df:   pandas df in standard format
        dist: (kb) distance around index SNP to clump
    Returns:
        pandas df with additional columns showing cluster number
    """
    df["cluster"] = None
    df = df.sort_values("pval", ascending=True)
    # Initiate
    clusnum = 1
    unclustered = pd.isnull(df["cluster"])
    # Continue clumping whilst there are unclustered SNPs
    while unclustered.any():
        # Get index row
        index_row = df.loc[unclustered, :].iloc[0]
        # Find other rows within set distance
        in_cluster = ( (df["chrom"] == index_row["chrom"]) &
                       (df["pos"] >= index_row["pos"] - dist * 1000) &
                       (df["pos"] <= index_row["pos"] + dist * 1000) &
                       unclustered )
        df.loc[in_cluster, "cluster"] = clusnum
        # Increase cluster number
        clusnum += 1
        unclustered = pd.isnull(df["cluster"])

    return df

def parse_args():
    """ Load command line args.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--inf', metavar="<file>", help=('Input file'), type=str, required=True)
    parser.add_argument('--out_common', metavar="<file>", help=('Output file for common vars >1% maf'), type=str, required=False)
    parser.add_argument('--out_rare', metavar="<file>", help=('Output file for rare vars <1% maf'), type=str, required=False)
    parser.add_argument('--out_combined', metavar="<file>", help=('Output file for all variants'), type=str, required=False)
    parser.add_argument('--window', metavar="<int>", help=('Distance window (kb) (default: 500)'), type=int, default=500)
    parser.add_argument('--pval', metavar="<float>", help=('P-value cutoff (default: 5e-8)'), type=float, default=5e-8)
    parser.add_argument('--exclMHC', help=('Exclude MHC region (GRCh37)'), action='store_true')
    parser.add_argument('--chrom', help=('Only do this chromosome'), type=str, required=False)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
