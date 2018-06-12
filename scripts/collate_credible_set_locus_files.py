#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import sys
import os
from glob import glob
import pandas as pd
import re
import argparse

def main():

    # Args
    in_pattern = "/nfs/users/nfs_e/em21/otcoregen/em21/finemapping/results/*/credible_set/*.cond.*.credible_sets.tsv"
    args = parse_args()

    # Load data
    study_files = glob(in_pattern)
    dfs = []
    for c, inf in enumerate(study_files):

        if c % 100 == 0:
            print("Loading {0} of {1}...".format(c + 1, len(study_files)))

        # Get trait and index snp name
        pattern = r"/nfs/users/nfs_e/em21/otcoregen/em21/finemapping/results/(.+)/credible_set/[0-9]+\.cond\.(.+)\.credible_sets\.tsv"
        mth_obj = re.match(pattern, inf)
        trait = mth_obj.group(1)
        index_snp = mth_obj.group(2)

        # Load study
        df = pd.read_csv(inf, sep="\t", header=0)
        df = df[df["is99_credset"] == 1]
        df.insert(0, "trait", trait)
        df.insert(1, "locus_index_snp", index_snp)
        dfs.append(df)

    # Merge studies into single df
    print("Merging datasets...")
    merged = pd.concat(dfs, axis=0, ignore_index=True)
    print(merged.shape)

    # Save
    print("Saving...")
    merged.to_csv(args.out, sep="\t", index=None, compression="gzip")

    return 0

def parse_args():
    """ Load command line args.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', metavar="<file>", help=('Output file'), type=str, required=False)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
