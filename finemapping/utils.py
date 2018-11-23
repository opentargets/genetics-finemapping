#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import dask.dataframe as dd
import numpy as np
import pandas as pd
from collections import OrderedDict

def load_sumstats(in_pq, study_id, cell_id=None, group_id=None, trait_id=None,
                  chrom=None, excl_mhc=None, min_maf=None, build='b37',
                  logger=None):
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

    # Create row-group filters
    row_grp_filters = [('study_id', '==', study_id)]
    if cell_id:
        row_grp_filters.append(('cell_id', '==', cell_id))
    if group_id:
        row_grp_filters.append(('group_id', '==', group_id))
    if trait_id:
        row_grp_filters.append(('trait_id', '==', trait_id))
    if chrom:
        row_grp_filters.append(('chrom', '==', str(chrom)))

    # Create column filters
    # TODO

    # Read file
    df = dd.read_parquet(in_pq,
                         filters=row_grp_filters,
                         engine='fastparquet')

    # Conversion to in-memory pandas
    df = df.compute(scheduler='single-threaded') # DEBUG
    df = df.astype(dtype=get_meta_info(type='sumstats'))

    # Apply row filters
    query = ' & '.join(
        ['{} {} "{}"'.format(filter[0], filter[1], filter[2])
         for filter in row_grp_filters] )
    df = df.query(query)

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

    # Extract EAF. TODO this needs to be changed to use estimated EAF if EAF_est
    if pd.isnull(df['eaf']).any():
        if logger:
            logger.warning('Warning: using MAF instead of EAF')
    df['eaf'] = np.where(pd.isnull(df['eaf']),
                             df['maf'],
                             df['eaf'])

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

def get_credset_out_columns():
    ''' Returns an OrderedDict to map credible set output columns.
        Useful for creating an empty credset df.
    Returns:
        OrderedDict
    '''
    return OrderedDict([
        ('study_id', 'study_id'),
        ('cell_id', 'cell_id'),
        ('gene_id', 'gene_id'),
        ('group_id', 'group_id'),
        ('trait_id', 'trait_id'),
        ('index_variant_id', 'variant_id_index'),
        ('variant_id', 'variant_id_tag'),
        ('pos', 'pos_tag'),
        ('chrom', 'chrom_tag'),
        ('ref_al', 'ref_al_tag'),
        ('alt_al', 'alt_al_tag'),
        ('beta', 'beta_tag'),
        ('se', 'se_tag'),
        ('pval', 'pval_tag'),
        ('beta_cond', 'beta_cond_tag'),
        ('se_cond', 'se_cond_tag'),
        ('pval_cond', 'pval_cond_tag'),
        ('logABF', 'logABF'),
        ('postprob', 'postprob'),
        ('postprob_cumsum', 'postprob_cumsum'),
        ('is95_credset', 'is95_credset'),
        ('is99_credset', 'is99_credset'),
        ('multisignal_method', 'multisignal_method')
    ])

def get_toploci_out_columns():
    ''' Returns an OrderedDict to map top loci set output columns.
        Useful for creating an empty credset df.
    Returns:
        OrderedDict
    '''
    return OrderedDict([
        ('study_id', 'study_id'),
        ('cell_id', 'cell_id'),
        ('gene_id', 'gene_id'),
        ('group_id', 'group_id'),
        ('trait_id', 'trait_id'),
        ('variant_id', 'variant_id'),
        ('chrom', 'chrom'),
        ('pos', 'pos'),
        ('ref_al', 'ref_al'),
        ('alt_al', 'alt_al'),
        ('beta', 'beta'),
        ('se', 'se'),
        ('pval', 'pval'),
        ('clump_method', 'clump_method')
    ])

def get_meta_info(type):
    ''' Returns a dict of meta data for dask
    Args:
        type [top_loci|cred_set|sumstats]
    '''
    if type == 'top_loci':
        meta = {
            'study_id': 'object',
            'cell_id': 'object',
            'gene_id': 'object',
            'group_id': 'object',
            'trait_id': 'object',
            'variant_id': 'object',
            # 'chrom': 'category',
            'chrom': 'object',
            'pos': 'int64',
            'ref_al': 'object',
            'alt_al': 'object',
            'beta': 'float64',
            'se': 'float64',
            'pval': 'float64',
            'clump_method': 'object'
        }
    elif type == 'cred_set':
        meta = {
            'study_id': 'object',
            'cell_id': 'object',
            'gene_id': 'object',
            'group_id': 'object',
            'trait_id': 'object',
            'variant_id_index': 'object',
            'variant_id_tag': 'object',
            'pos_tag': 'int64',
            # 'chrom_tag': 'category',
            'chrom_tag': 'object',
            'ref_al_tag': 'object',
            'alt_al_tag': 'object',
            'beta_tag': 'float64',
            'se_tag': 'float64',
            'pval_tag': 'float64',
            'beta_cond_tag': 'float64',
            'se_cond_tag': 'float64',
            'pval_cond_tag': 'float64',
            'logABF': 'float64',
            'postprob': 'float64',
            'postprob_cumsum': 'float64',
            'is95_credset': 'bool',
            'is99_credset': 'bool',
            'multisignal_method': 'object'
        }
    elif type == 'sumstats':
        meta = {
            'variant_id_b37': 'object',
            'pos_b37': 'int64',
            'ref_al': 'object',
            'alt_al': 'object',
            'beta': 'float64',
            'se': 'float64',
            'pval': 'float64',
            'n_samples_variant_level': 'float64',
            'n_samples_study_level': 'float64',
            'n_cases_variant_level': 'float64',
            'n_cases_study_level': 'float64',
            'eaf': 'float64',
            'maf': 'float64',
            'info': 'float64',
            'is_cc': 'bool',
            'study_id': 'object',
            'group_id': 'object',
            'cell_id': 'object',
            'gene_id': 'object',
            'trait_id': 'object',
            'chrom': 'object'
        }
    return meta

def is_local(path):
    ''' Checks if a path is local '''
    if '://' in path:
        return False
    else:
        return True

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
