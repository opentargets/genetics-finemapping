#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import dask.dataframe as dd
import pandas as pd
from collections import OrderedDict
import re


def load_sumstats(in_pq, study_id, phenotype_id=None, bio_feature=None,
                  chrom=None, excl_mhc=None, min_maf=None, logger=None):
    ''' Loads summary statistics from Open Targets parquet format using dask:
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
    if bio_feature:
        row_grp_filters.append(('bio_feature', '==', bio_feature))
    if phenotype_id:
        row_grp_filters.append(('phenotype_id', '==', phenotype_id))
    if chrom:
        # HACK: This workaround is needed because for molecular trait studies we
        # partition the files by chromosome. But spark/dask automatically infer the
        # column type as integer, whereas for GWAS this column is in the metadata
        # as a string. So there is no way to give a single column type in this row
        # group filter.
        if re.search('molecular_trait', in_pq):
            row_grp_filters.append(('chrom', '==', int(chrom)))
        else:
            row_grp_filters.append(('chrom', '==', str(chrom)))
        # Note: without this fix, you get errors like:
        # pyarrow.lib.ArrowNotImplementedError: Function equal has no kernel matching input types (array[int32], scalar[string])

    # Create column filters
    cols_to_keep = ['study_id', 'phenotype_id', 'bio_feature', 'chrom', 'pos',
                    'ref', 'alt', 'beta', 'se', 'pval', 'n_total', 'n_cases',
                    'eaf', 'is_cc']

    # Read file
    # gather_statistics=True seems to be necessary, otherwise get a dask
    # error "Cannot apply filters with gather_statistics=False"
    # index=False is also essential to avoid dask errors with newer dask versions
    df = dd.read_parquet(in_pq,
                         columns=cols_to_keep,
                         filters=row_grp_filters,
                         gather_statistics=True,
                         engine='pyarrow',
                         index=False)

    # Conversion to in-memory pandas
    df = (
        df.compute(scheduler='single-threaded')
          .astype(dtype=get_meta_info(type='sumstats'))
    )

    # Apply row filters
    if study_id:
        df = df.loc[df['study_id'] == study_id, :]
    if bio_feature:
        df = df.loc[df['bio_feature'] == bio_feature, :]
    if phenotype_id:
        df = df.loc[df['phenotype_id'] == phenotype_id, :]
    if chrom:
        df = df.loc[df['chrom'] == chrom, :]

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
                       (df['pos'] >= 28477797) &
                       (df['pos'] <= 33448354) )
            df = df.loc[~is_mhc, :]
        elif excl_mhc == 'b38':
            is_mhc = ( (df['chrom'] == '6') &
                       (df['pos'] >= 28510120) &
                       (df['pos'] <= 33480577) )
            df = df.loc[~is_mhc, :]
    
    # Create a variant ID
    df['variant_id'] = (
        df.loc[:, ['chrom', 'pos', 'ref', 'alt']]
        .apply(lambda row: ':'.join([str(x) for x in row]), axis=1)
    )

    return df

# import pyspark.sql
# from pyspark.sql.functions import col
# def load_sumstats(in_pq, study_id, phenotype_id=None, bio_feature=None,
#                   chrom=None, excl_mhc=None, min_maf=None, logger=None):
#     ''' Loads summary statistics from Open Targets parquet format using Spark
#     '''

#     # Get spark session
#     spark = (
#         pyspark.sql
#         .SparkSession
#         .builder
#         .master("local")
#         .getOrCreate()
#     )

#     # Create column filters
#     cols_to_keep = ['study_id', 'phenotype_id', 'bio_feature', 'chrom', 'pos',
#                     'ref', 'alt', 'beta', 'se', 'pval', 'n_total', 'n_cases',
#                     'eaf', 'is_cc']

#     # Load
#     df_spark = (
#         spark.read.parquet(in_pq)
#         .select(cols_to_keep)
#     )

#     # Apply required filters
#     df_spark = (
#         df_spark.filter(
#             (col('study_id') == study_id) &
#             (col('chrom') == chrom)
#         )
#     )

#     # Apply optional filters
#     if phenotype_id:
#         df_spark = df_spark.filter(col('phenotype_id') == phenotype_id)
#     if bio_feature:
#         df_spark = df_spark.filter(col('bio_feature') == bio_feature)

#     # Convert to pandas
#     df = df_spark.toPandas()

#     # Exclude on MAF
#     if min_maf:
#         to_exclude = (df['eaf'].apply(eaf_to_maf) < min_maf)
#         df = df.loc[~to_exclude, :]
    
#     # Exclude MHC
#     if excl_mhc and (df.chrom == '6').any():
#         # Exclude MHC
#         if excl_mhc == 'b37':
#             is_mhc = ( (df['chrom'] == '6') &
#                        (df['pos'] >= 28477797) &
#                        (df['pos'] <= 33448354) )
#             df = df.loc[~is_mhc, :]
#         elif excl_mhc == 'b38':
#             is_mhc = ( (df['chrom'] == '6') &
#                        (df['pos'] >= 28510120) &
#                        (df['pos'] <= 33480577) )
#             df = df.loc[~is_mhc, :]

#     # Create a variant ID
#     df['variant_id'] = (
#         df.loc[:, ['chrom', 'pos', 'ref', 'alt']]
#         .apply(lambda row: ':'.join([str(x) for x in row]), axis=1)
#     )

#     return df

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
        ('bio_feature', 'bio_feature'),
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
        ('bio_feature', 'bio_feature'),
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
            'bio_feature': 'object',
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
            'bio_feature': 'object',
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
            'bio_feature': 'object',
            'chrom': 'str',
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
    elif type == 'finemap_snp':
        meta = {
            'locus_name': 'str',
            'index': 'int64',
            'rsid': 'str',
            'chromosome': 'str',
            'position': 'int64',
            'allele1': 'str',
            'allele2': 'str',
            'maf': 'float64',
            'beta': 'float64',
            'se': 'float64',
            'z': 'float64',
            'prob': 'float64',
            'log10bf': 'float64',
            'mean': 'float64',
            'sd': 'float64',
            'mean_incl': 'float64',
            'sd_incl': 'float64'
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

def df_empty(columns, dtypes, index=None):
    ''' Creat an empty df
    '''
    assert len(columns)==len(dtypes)
    df = pd.DataFrame(index=index)
    for c, d in zip(columns, dtypes):
        df[c] = pd.Series(dtype=d)
    return df
