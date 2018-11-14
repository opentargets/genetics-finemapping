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

def main():

    # Args
    row_offset = 500000 # Approx num rows per row-group

    #
    # Process an example Neale study
    #

    print('Making example Neale...')

    # Download Neale 50 (height)
    url = 'gs://genetics-portal-sumstats/gwas/genome_wide/NEALEUKB_50/UKB_50/*-NEALEUKB_50-UKB_50.tsv.gz'
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

    #
    # Process an example GTEX tissue
    #

    print('Making example GTEX...')

    # Download GTEX whole blood (UBERON_0000178)
    url = 'gs://genetics-portal-sumstats/molecular_qtl/eqtl/GTEX7/UBERON_0000178/*/*.tsv.gz'
    dfs = []
    for inf in glob_gcs(url, head=500):
        # Get information about file
        _, study_id, cell_id, gene_id = os.path.split(inf)[-1].replace('.tsv.gz', '').split('-')
        # Load dataset
        df = dd.read_table(inf, sep='\t', compression='gzip', blocksize=None)
        df['study_id'] = study_id
        df['cell_id'] = cell_id
        df['gene_id'] = gene_id
        df['group_id'] = gene_id
        df['trait_id'] = 'eqtl'
        dfs.append(df.compute())

    # Merge
    merged = pd.concat(dfs)
    merged = merged.astype(dtype={'chrom':'str'})
    print(merged.info())

    # Sort
    merged = merged.sort_values(
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
    os.makedirs('molecular_qtl', exist_ok=True)
    out_path = 'molecular_qtl/GTEX7'
    merged.to_parquet(
        out_path,
        engine='fastparquet',
        compression='snappy',
        file_scheme='hive',
        row_group_offsets=row_offset
    )

    # Save example as tsv
    out_path = 'molecular_qtl/GTEX7.tsv.gz'
    merged.head(1000).to_csv(
        out_path,
        sep='\t',
        index=None,
        compression='gzip'
    )

    return 0

def glob_gcs(url, head=None):
    ''' Uses gsutil to glob a gcs pattern
    '''
    # Make the command
    cmd = 'gsutil ls {0}'.format(url)
    if head:
        cmd += ' | head -{0}'.format(head)
    # Get result
    p = sp.Popen(cmd, stdout=sp.PIPE, shell=True)
    result = p.communicate()[0].decode().split('\n')
    return(result[:-1])


if __name__ == '__main__':

    main()
