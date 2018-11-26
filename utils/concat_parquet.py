#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import argparse
import pandas as pd

def main():

    # Args
    args = parse_args()

    # Load
    dfs = (pd.read_parquet(inf, engine='fastparquet')
           for inf in args.in_parquets)

    # Concatenate
    full_df = pd.concat(dfs, ignore_index=True)

    # Write
    full_df.to_parquet(
        args.out,
        engine='fastparquet',
        compression='snappy',
        row_group_offsets=500000
    )

    return 0

def parse_args():
    ''' Load command line args
    '''
    p = argparse.ArgumentParser()

    # Add input files
    p.add_argument('--in_parquets',
                   metavar="<file>",
                   help=("List of parquet files to concatenate"),
                   type=str,
                   nargs='+',
                   required=True)
    p.add_argument('--out',
                   metavar="<file>",
                   help=("Concatenated parquet file"),
                   type=str,
                   required=True)

    args = p.parse_args()

    return args

if __name__ == '__main__':

    main()
