#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import dask.dataframe as dd
import numpy as np
import pandas as pd
from collections import OrderedDict
import os

def load_sumstats(in_pq, study_id, phenotype_id=None, biofeature=None,
                  chrom=None, excl_mhc=None, min_maf=None, logger=None):
    ''' Loads summary statistics from Open Targets parquet format:
        - Loads only required rows
        - Converts to pandas
        - TODO only extract required columns
        - TODO extract eaf
    Args:
        excl_mhc (b37|b38|None): whether to exclude the MHC region
    '''

    #
    # Load
    #

    # Create row-group filters
    row_grp_filters = [('study_id', '==', study_id)]
    if phenotype_id:
        row_grp_filters.append(('phenotype_id', '==', phenotype_id))
    if chrom:
        row_grp_filters.append(('chrom', '==', str(chrom)))

    # Add biofeature to path
    if biofeature:
        in_pq = os.path.join(in_pq, 'biofeature={}'.format(biofeature))

    # Create column filters
    cols_to_keep = ['study_id', 'phenotype_id', 'chrom', 'pos',
                    'ref', 'alt', 'beta', 'se', 'pval', 'n_total', 'n_cases',
                    'eaf', 'is_cc']

    # Read file
    in_pq_pattern = os.path.join(in_pq, '*.parquet')
    df = dd.read_parquet(in_pq_pattern,
                         columns=cols_to_keep,
                         filters=row_grp_filters,
                         engine='fastparquet')

    # Conversion to in-memory pandas
    df = df.compute(scheduler='single-threaded')
    df = df.astype(dtype=get_meta_info(type='sumstats'))

    # Apply row filters
    query_parts = []
    for part in row_grp_filters:
        if isinstance(part[2], int):
            query_parts.append('{} {} {}'.format(part[0], part[1], part[2]))
        else:
            query_parts.append('{} {} "{}"'.format(part[0], part[1], part[2]))
    query = ' & '.join(query_parts)
    df = df.query(query)

    # Add biofeature back in
    df.loc[:, 'biofeature'] = biofeature

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
    
    # Create a variant ID
    df['variant_id'] = (
        df.loc[:, ['chrom', 'pos', 'ref', 'alt']]
        .apply(lambda row: ':'.join([str(x) for x in row]), axis=1)
    )

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
        ('phenotype_id', 'phenotype_id'),
        ('biofeature', 'biofeature'),
        ('lead_variant_id', 'lead_variant_id'),
        ('lead_chrom', 'lead_chrom'),
        ('lead_pos', 'lead_pos'),
        ('lead_ref', 'lead_ref'),
        ('lead_alt', 'lead_alt'),
        ('variant_id', 'tag_variant_id'),
        ('chrom', 'tag_chrom'),
        ('pos', 'tag_pos'),
        ('ref', 'tag_ref'),
        ('alt', 'tag_alt'),
        ('beta', 'tag_beta'),
        ('se', 'tag_se'),
        ('pval', 'tag_pval'),
        ('beta_cond', 'tag_beta_cond'),
        ('se_cond', 'tag_se_cond'),
        ('pval_cond', 'tag_pval_cond'),
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
        ('phenotype_id', 'phenotype_id'),
        ('biofeature', 'biofeature'),
        ('variant_id', 'variant_id'),
        ('chrom', 'chrom'),
        ('pos', 'pos'),
        ('ref', 'ref'),
        ('alt', 'alt'),
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
            'phenotype_id': 'object',
            'biofeature': 'object',
            'variant_id': 'object',
            'chrom': 'object',
            'pos': 'int64',
            'ref': 'object',
            'alt': 'object',
            'beta': 'float64',
            'se': 'float64',
            'pval': 'float64',
            'clump_method': 'object'
        }
    elif type == 'cred_set':
        meta = {
            'study_id': 'object',
            'phenotype_id': 'object',
            'biofeature': 'object',
            'lead_variant_id': 'object',
            'lead_chrom': 'object',
            'lead_pos': 'int64',
            'lead_ref': 'object',
            'lead_alt': 'object',
            'tag_variant_id': 'object',
            'tag_chrom': 'object',
            'tag_pos': 'int64',
            'tag_ref': 'object',
            'tag_alt': 'object',
            'tag_beta': 'float64',
            'tag_se': 'float64',
            'tag_pval': 'float64',
            'tag_beta_cond': 'float64',
            'tag_se_cond': 'float64',
            'tag_pval_cond': 'float64',
            'logABF': 'float64',
            'postprob': 'float64',
            'postprob_cumsum': 'float64',
            'is95_credset': 'bool',
            'is99_credset': 'bool',
            'multisignal_method': 'object'
        }
    elif type == 'sumstats':
        meta = {
            'study_id': 'object',
            'phenotype_id': 'object',
            # 'biofeature': 'object',
            'chrom': 'object',
            'pos': 'int64',
            'ref': 'object',
            'alt': 'object',
            'beta': 'float64',
            'se': 'float64',
            'pval': 'float64',
            'n_total': 'float64',
            'n_cases': 'float64',
            'eaf': 'float64',
            'is_cc': 'bool'
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
