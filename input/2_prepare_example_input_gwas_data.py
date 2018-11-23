#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import sys
import os
import pandas as pd
from dask import dataframe as dd
from dask import delayed
import subprocess as sp
import dask

def main():

    #
    # Process an example Neale study
    #

    row_offset = 500000 # Approx num rows per row-group

    print('Making example Neale...')

    # Download Neale 50 (height)
    url = 'raw/NEALEUKB_50/UKB_50/*-NEALEUKB_50-UKB_50.tsv.gz'
    df = dd.read_table(url, sep='\t', compression='gzip', blocksize=None).compute()
    df['study_id'] = 'NEALEUKB_50'
    df['cell_id'] = ''
    df['gene_id'] = ''
    df['group_id'] = ''
    df['trait_id'] = 'UKB_50'

    # Set types
    df = df.astype(dtype={'chrom':'str'})
    print(df.info())

    # Sort
    df = df.sort_values(
        [
            'study_id',
            'trait_id',
            'cell_id',
            'group_id',
            'chrom',
            'pos_b37'
        ]
    )

    # Save as parquet
    os.makedirs('gwas', exist_ok=True)
    out_path = 'gwas/NEALEUKB_50'
    df.to_parquet(
        out_path,
        engine='fastparquet',
        compression='snappy',
        file_scheme='hive',
        row_group_offsets=row_offset
    )

    # Save example as tsv
    out_path = 'gwas/NEALEUKB_50.tsv.gz'
    df.head(1000).to_csv(
        out_path,
        sep='\t',
        index=None,
        compression='gzip'
    )

    return 0

if __name__ == '__main__':

    main()
