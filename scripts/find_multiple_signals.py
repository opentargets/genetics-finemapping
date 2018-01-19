#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#
"""
Takes *.jma output from GCTA cojo

Usage: python find_multiple_signals.py <inf> <kb> <outf>

Args:
    inf:  input jma file
    kb:   distance on which to cluster SNPs at the same locus (e.g. 500kb)
    outf: output file

Output, 2 columns:
    index_snp: index snp for signal
    conditional_snps: other snps that it should be conditional on in col2

"""

import sys
import pandas as pd

def main():

    # Args
    inf = sys.argv[1]
    cluster_kb = int(sys.argv[2])
    outfile = sys.argv[3]

    # Load data to pandas
    data = pd.read_csv(inf, header=0, sep="\t")
    # print(data.head())

    # Find loci with multiple signals
    multi_loci_dict = make_multi_loci_dict(data, cluster_kb)

    # Write conditional config file
    with open(outfile, "w") as out_h:
        # Write header
        out_h.write("\t".join(["index_snp", "conditional_snps"]) + "\n")
        # Write rows
        for key, value in multi_loci_dict.items():
            out_h.write("\t".join([key, ";".join(value)]) + "\n")

    return 0

def make_multi_loci_dict(df, dist=500):
    """ For each SNP, find other variants on the same chromosome that are within
        a set distance from it.
    Args:
        df:   pandas df of GCTA jma data
        dist: distance in Kb to be considered a shared locus
    Returns:
        dict {index SNP: [shared locus snps, ...]}
    """
    d = {}
    for _, index_row in df.iterrows():
        is_shared = ( (df["Chr"] == index_row["Chr"]) &
                      (df["SNP"] != index_row["SNP"]) &
                      ((df["bp"] - index_row["bp"]).apply(abs) <= dist * 1000) )
        d[index_row["SNP"]] = df.loc[is_shared, "SNP"].tolist()
    return d

if __name__ == '__main__':

    main()
