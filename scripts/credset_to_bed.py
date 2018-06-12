#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#
"""
Reads the output of multiple credible set output files and merges them together
into a BED file (0-indexed).

Output cols:
    chrom
    start bp (0-indexed)
    end bp (not inclusive)
    name: index snp
    size: number of snps in cred set
    Credible set SNPs
"""

import sys
import pandas as pd
import argparse
import re

def main():

    # Args
    args = parse_args()
    cred_thresh = float(args.credset) / 100

    # Initiate out data
    out_rows = []

    # Process each locus
    for inf in args.infiles:

        # Get index snp
        mth_obj = re.search("\.cond\.(rs[0-9]+)\.credible_sets\.", inf)
        index_snp = mth_obj.group(1)

        # Load and sort
        data = pd.read_csv(inf, sep="\t", header=0).sort_values("postprob_cumsum")

        # Find rows in credible set
        set_idx = data["postprob_cumsum"].gt(cred_thresh).tolist().index(True)
        data_cred = data.iloc[0:set_idx + 1, :]

        # Make out row
        out_row = [data_cred.iloc[0, 0],
                   min(data_cred.bp) - 1,
                   max(data_cred.bp) + 1,
                   index_snp,
                   data_cred.shape[0],
                   ";".join(data_cred.SNP.tolist())]
        out_rows.append(out_row)

    # Make a df
    out_df = pd.DataFrame(out_rows,
                          columns=["chrom", "start", "end", "index", "size", "set"])
    out_df = out_df.sort_values(["chrom", "start", "end"])
    out_df.to_csv(args.outbed, sep="\t", index=None, header=None)

    return 0

def parse_args():
    """ Load command line args.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--infiles', metavar="<file1> <file2> ...", help=('Input files'), nargs="+", type=str, required=True)
    parser.add_argument('--outbed', metavar="<file>", help=("Output bed file"), type=str, required=True)
    parser.add_argument('--credset', metavar="<int>", help=('Credible set % threshold (default: 95)'), default=95, type=int, required=False)

    args = parser.parse_args()

    return args

if __name__ == '__main__':

    main()
