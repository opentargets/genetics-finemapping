#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import dask.dataframe as dd
import numpy as np
import pandas as pd

def load_sumstats(in_pq, study_id, cell_id=None, gene_id=None,
                  group_id=None, chrom=None, excl_mhc=None, min_maf=None,
                  build='b37'):
    ''' Loads summary statistics from Open Targets parquet format:
        - Loads only required rows
        - Converts to pandas
        - Extracts N samples, N cases, EAF
        - TODO only extract required columns
        - TODO extract eaf
    Args:
        excl_mhc (b37|b38|None): whether to exclude the MHC region
        build (b37|b38): which build to use
    '''

    #
    # Load
    #

    # Create row filters
    row_filters = [('study_id', '==', study_id)]
    if cell_id:
        row_filters.append(('cell_id', '==', cell_id))
    if gene_id:
        row_filters.append(('gene_id', '==', gene_id))
    if group_id:
        row_filters.append(('group_id', '==', group_id))
    if chrom:
        row_filters.append(('chrom', '==', str(chrom)))

    # Create column filters
    # TODO

    # Read file
    df = dd.read_parquet(in_pq,
                         filters=row_filters,
                         engine='fastparquet')

    # Conversion to in-memory pandas
    df = df.compute()

    #
    # Extract fields
    #

    # Extract n_samples
    df['n_samples'] = np.where(pd.isnull(df['n_samples_variant_level']),
                               df['n_samples_study_level'],
                               df['n_samples_variant_level'])
    # Extract n_cases
    df['n_cases'] = np.where(pd.isnull(df['n_cases_variant_level']),
                             df['n_cases_study_level'],
                             df['n_cases_variant_level'])

    # Extract required build
    df['variant_id'] = df['variant_id_{0}'.format(build)]
    df['pos'] = df['pos_{0}'.format(build)]

    #
    # Make exclusions
    #

    # Exclude on MAF
    if min_maf:
        to_exclude = ( df['eaf'].apply(eaf_to_maf) < min_maf )
        df = df.loc[~to_exclude, :]

    # Exclude MHC
    if excl_mhc and (df.chrom == '6').any():
        # Exclude MHC
        if excl_mhc == 'b37':
            is_mhc = ( (df['chrom'] == '6') &
                       (df['pos_b37'] >= 28477797) &
                       (df['pos_b37'] <= 33448354) )
            df = df.loc[~is_mhc, :]
        elif excl_mhc == 'b38':
            is_mhc = ( (df['chrom'] == '6') &
                       (df['pos_b38'] >= 28510120) &
                       (df['pos_b38'] <= 33480577) )
            df = df.loc[~is_mhc, :]
    return df

def eaf_to_maf(eaf):
    ''' Convert effect allele frequency to MAF
    '''
    eaf = float(eaf)
    if eaf <= 0.5:
        return eaf
    else:
        return 1 - eaf

def extract_window(sumstats, chrom, pos, window):
    ''' Extracts a window around a genomic position
    Args:
        sumstats (pd.df)
        chrom (str)
        pos (int)
        window (int): kb around pos to extract
    Returns:
        pd.df
    '''
    in_window = ( (sumstats['chrom'] == chrom) &
                  (sumstats['pos'] >= pos - 1000 * window) &
                  (sumstats['pos'] <= pos + 1000 * window) )
    sumstat_wind = sumstats.loc[in_window, :]

    return sumstat_wind


# def load_manifest(in_data):
#     ''' Loads a dataframe containing information about partitions on which to
#         run the finemapping pipeline.
#     Args:
#         in_data (parquet): parquet file containing sumstats
#     Returns:
#         dd
#     '''
#     req_cols = ['study_id', 'cell_id', 'gene_id', 'trait_id', 'chrom']
#     df = ( dd.read_parquet(in_data,
#                            columns=req_cols,
#                            filters=[('chrom', 'in', ['22'])],
#                            # filters=[('chrom', 'in', ['21', '22'])],
#                            engine='fastparquet')
#              .drop_duplicates() )
#     return df
