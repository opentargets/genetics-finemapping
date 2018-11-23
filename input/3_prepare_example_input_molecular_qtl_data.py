#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import os
import sys
from dask import dataframe as dd
from collections import OrderedDict
from dask.diagnostics import ProgressBar

def main():

    #
    # Process an example GTEX tissue
    #

    # Download GTEX whole blood (UBERON_0000178)
    inf = 'raw/GTEX/UBERON_0000178/Whole_Blood.allpairs.txt'

    # Load
    df = dd.read_table(inf, sep='\t')

    # Extract variant information
    df.variant_id = df.variant_id.str.replace('_b37', '')

    # Split chrom_pos_ref_alt
    df = df.assign(
        chrom_pos_ref_alt = df.variant_id.map(lambda x: x.split('_'))
    )

    # Assign to different columns
    print('Assigning variant columns...')
    df = df.assign(
        chrom=df.chrom_pos_ref_alt.map(lambda x: x[0], meta=('chrom', 'str')),
        pos=df.chrom_pos_ref_alt.map(lambda x: x[1], meta=('pos', 'int')),
        ref=df.chrom_pos_ref_alt.map(lambda x: x[2], meta=('ref', 'str')),
        alt=df.chrom_pos_ref_alt.map(lambda x: x[3], meta=('alt', 'str'))
    )

    # Clean gene name
    df = df.assign(
        gene_clean = df.gene_id.map(lambda x: x.split('.')[0],
                                    meta=('gene_clean', 'str'))
        )

    # Add sample size
    df['n_samples_study_level'] = 369

    # Set whether study is case-control
    df['is_cc'] = False

    # Add other needed columns
    df['study_id'] = 'GTEX7'
    df['cell_id'] = 'UBERON_0000178'
    df['trait_id'] = 'eqtl'
    df['gene_id'] = df['gene_clean']
    df['group_id'] = df['gene_clean']

    # Rename columns
    req_cols = OrderedDict([
        ('variant_id', 'variant_id_b37'),
        ('chrom', 'chrom'),
        ('pos', 'pos_b37'),
        ('ref', 'ref_al'),
        ('alt', 'alt_al'),
        ('slope', 'beta'),
        ('slope_se', 'se'),
        ('pval_nominal', 'pval'),
        ('n', 'n_samples_variant_level'),
        ('n_samples_study_level', 'n_samples_study_level'),
        ('n_cases_variant_level', 'n_cases_variant_level'),
        ('ncases', 'n_cases_study_level'),
        ('eaf', 'eaf'),
        ('maf', 'maf'),
        ('info', 'info'),
        ('is_cc', 'is_cc'),
        ('study_id', 'study_id'),
        ('cell_id', 'cell_id'),
        ('trait_id', 'trait_id'),
        ('gene_id', 'gene_id'),
        ('group_id', 'group_id')
    ])
    df = df.rename(columns=req_cols)
    df = df.loc[:, req_cols.values()]

    # Set dtypes
    dtypes_dict = {
        'variant_id_b37': str,
        'chrom': str,
        'pos_b37': int,
        'ref_al': str,
        'alt_al': str,
        'beta': float,
        'se': float,
        'pval': float,
        # 'n_samples_variant_level': int,
        'n_samples_study_level': int,
        # 'n_cases_variant_level': int,
        # 'n_cases_study_level': int,
        'eaf': float,
        'maf': float,
        'info': float,
        'is_cc': bool,
        'gene_id': str,
        'study_id': str,
        'cell_id': str,
        'trait_id': str,
        'group_id': str
    }
    df = df.astype(dtype=dtypes_dict)

    # Save as parquet
    print('Saving...')
    os.makedirs('molecular_qtl', exist_ok=True)
    out_path = 'molecular_qtl/GTEX7_fastparquet_partitioned'
    with ProgressBar():
        df.to_parquet(
            out_path,
            # engine='pyarrow',
            engine='fastparquet',
            compression='snappy',
            partition_on=['group_id']
            # file_scheme='hive'
            # row_group_offsets=row_offset
        )

    # Save example as tsv
    out_path = 'molecular_qtl/GTEX7.tsv.gz'
    df.head(1000).to_csv(
        out_path,
        sep='\t',
        index=None,
        compression='gzip'
    )


    # dfs = []
    # print('Getting list of files...')
    # in_files = glob_gcs(url)
    # # in_files = glob_gcs(url, head=500)
    # print('Making df list...')
    # for i, inf in enumerate(in_files):
    #     print('Adding file {0}...'.format(i))
    #     # Get information about file
    #     _, study_id, cell_id, gene_id = os.path.split(inf)[-1].replace('.tsv.gz', '').split('-')
    #     # Load dataset
    #     df = dd.read_table(inf, sep='\t', compression='gzip', blocksize=None)
    #     df['study_id'] = study_id
    #     df['cell_id'] = cell_id
    #     df['gene_id'] = gene_id
    #     df['group_id'] = gene_id
    #     df['trait_id'] = 'eqtl'
    #     df['chrom'] = df['chrom'].astype(str)
    #
    #     dfs.append(df)
    #
    # # Merge
    # print('Concatenating dfs...')
    # merged = dd.concat(dfs)
    #
    # # Repartition
    # print('Repartitioning dfs...')
    # nparts = min(100, merged.npartitions)
    # merged = merged.repartition(npartitions=nparts)
    #
    # # Sort
    # # print('Sorting...')
    # # merged = merged.sort_values(
    # #     [
    # #         'study_id',
    # #         'trait_id',
    # #         'cell_id',
    # #         'group_id',
    # #         'chrom',
    # #         'pos_b37'
    # #     ]
    # # )
    #
    # # Save as parquet
    # print('Saving...')
    # os.makedirs('molecular_qtl', exist_ok=True)
    # out_path = 'molecular_qtl/GTEX7'
    # merged.to_parquet(
    #     out_path,
    #     engine='fastparquet',
    #     compression='snappy'
    #     # file_scheme='hive'
    #     # row_group_offsets=row_offset
    # )
    #
    # # Save example as tsv
    # out_path = 'molecular_qtl/GTEX7.tsv.gz'
    # merged.head(1000).to_csv(
    #     out_path,
    #     sep='\t',
    #     index=None,
    #     compression='gzip'
    # )

    return 0



if __name__ == '__main__':

    main()
