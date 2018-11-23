#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import gcta as fm_gcta
import utils as fm_utils
import os
import pandas as pd
from collections import OrderedDict

def detect_top_loci(sumstats, in_plink, temp_dir,
        method='conditional',
        maf=0.01,
        cojo_p=5e-8,
        cojo_window=500,
        cojo_collinear=0.9,
        clump_dist=500,
        clump_p=5e-8,
        logger=None):
    ''' Find top loci for a given study
    Args:
        TODO
        method ([conditional|distance]): method for identifying top loci
    '''

    temp_dir = os.path.join(temp_dir, 'top_loci')

    # Calculate how many variants below p threshold there are
    rows_below_threshold = (sumstats['pval'] <= cojo_p).sum()
    if logger:
        logger.info('{0} variants with p < {1}'.format(
            rows_below_threshold, cojo_p
        ))

    # Skip if clumping there are no variants with pval < p_threshold
    if rows_below_threshold == 0:
        if logger:
            logger.info('No variants with pval < p_threshold, skipping clumping...')
        # Use empty sumstats df
        top_loci = sumstats.head(0)
    # Detect top loci using GCTA-cojo
    elif method == 'conditional':
        if logger:
            logger.info('Clumping using GCTA-cojo...')
        # Use GCTA-cojo to perform conditional analysis
        top_loci = fm_gcta.get_conditional_top_loci(
            sumstats,
            in_plink,
            temp_dir,
            maf=maf,
            cojo_p=cojo_p,
            cojo_window=cojo_window,
            cojo_collinear=cojo_collinear,
            logger=logger
        )
    # Cluster using distance based clumping
    elif method == 'distance':
        if logger:
            logger.info('Clumping by distance...')
        top_loci = get_distance_top_loci(
            sumstats,
            clump_dist=clump_dist,
            clump_p=clump_p
        )
    # Raise error if conditional or distance is not selected
    else:
        if logger:
            logger.error('Method must be one of [conditional|distance]')
        raise ArgumentError('Method must be one of [conditional|distance]')

    # Add column specifying method used
    top_loci['clump_method'] = method

    # Format output table
    top_loci = format_top_loci_output(top_loci)

    return top_loci

def format_top_loci_output(top_loci):
    ''' Formats the top_loci table for output
    Args:
        top_loci (pd.df)
    Returns:
        pd.df
    '''
    cols = fm_utils.get_toploci_out_columns()
    meta = fm_utils.get_meta_info(type='top_loci')
    return top_loci.loc[:, cols.keys()] \
                   .rename(columns=cols) \
                   .astype(dtype=meta)

def get_distance_top_loci(sumstats, clump_dist=500, clump_p=5e-8):
    ''' Clump top loci based on distance.
    Args:
        sumstats (pd.df): summary stats
        clump_dist (int): assume loci within this distance (kb) are same signal
        clump_p (float): max p-value for a top snp
    Returns:
        sumstat table containing only top loci (pd.df)
    '''

    # Filter by p-value
    ss = sumstats.loc[sumstats['pval'] <= clump_p, :].copy()

    # Perform clustering
    ss.loc[:, 'cluster'] = None
    # ss.loc[:, "cluster"] = None
    ss["cluster"] = None
    ss = ss.sort_values("pval", ascending=True)
    # Initiate
    clusnum = 1
    unclustered = pd.isnull(ss["cluster"])
    # Continue clumping whilst there are unclustered SNPs
    while unclustered.any():
        # Get index row
        index_row = ss.loc[unclustered, :].iloc[0]
        # Find other rows within set distance
        in_cluster = (
            (ss["chrom"] == index_row["chrom"]) &
            (ss["pos"] >= index_row["pos"] - clump_dist * 1000) &
            (ss["pos"] <= index_row["pos"] + clump_dist * 1000) &
            unclustered
        )
        ss.loc[in_cluster, "cluster"] = clusnum
        # Increase cluster number
        clusnum += 1
        unclustered = pd.isnull(ss["cluster"])

    # Extract only top loci
    selected_snps = ss.drop_duplicates(subset='cluster', keep='first') \
                      .variant_id_b37
    top_loci = sumstats.loc[sumstats['variant_id'].isin(selected_snps), :]

    return top_loci


class ArgumentError(Exception):
    ''' Error to raise if invalid argument is used
    '''
    pass
