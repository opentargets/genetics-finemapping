#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#
"""
Reads the credible set bed files provided and outputs simple stat summary and
histogram.
"""

import sys
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def main():

    # Args
    args = parse_args()
    names = args.names if args.names else list(range(1, len(args.infiles) + 1))

    # Load all datasets into a single df
    dfs = []
    for idx, inf in enumerate(args.infiles):
        df = pd.read_csv(inf, sep="\t", header=None)
        df.columns = ["chrom", "start", "end", "index", "size", "set"]
        df["Group"] = names[idx]
        dfs.append(df)
    data = pd.concat(dfs)

    # Output mean and median
    data_stats = data.groupby('Group').agg({'size': [np.mean, np.median]})
    data_stats.to_csv(args.out, sep="\t")

    # # Make plot
    # if args.plot:
    #
    #     g = sns.FacetGrid(data, row="Group", margin_titles=True)
    #     g = g.map(plt.hist, "size", bins=100, range=(0, 100))
    #     g.savefig(args.plot)






    return 0

def parse_args():
    """ Load command line args.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--infiles', metavar="<file1> <file2> ...", help=('Input files'), nargs="+", type=str, required=True)
    parser.add_argument('--out', metavar="<file>", help=("Output file for stats"), type=str, required=True)
    parser.add_argument('--plot', metavar="<file>", help=("Output file for histogram"), type=str, required=False)
    parser.add_argument('--names', metavar="<str> <str> ...", help=("Names for datasets"), nargs="+", type=str, required=False)

    args = parser.parse_args()

    return args

if __name__ == '__main__':

    main()
