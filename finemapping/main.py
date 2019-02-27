#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import utils as fm_utils
import top_loci as fm_top_loci
import credible_set as fm_credible_set
import tempfile
import pandas as pd

def run_single_study(in_pq,
                     in_plink,
                     study_id,
                     phenotype_id=None,
                     biofeature=None,
                     chrom=None,
                     analysis_config=None,
                     tmp_dir=tempfile.gettempdir(),
                     method='conditional',
                     pval_threshold=5e-8,
                     logger=None):
    ''' Runs the top loci and credible set analysis on a single study
    '''

    pd.options.mode.chained_assignment = None # Silence pandas warning

    # TEMP solution (it is reqired that chrom is not None)
    assert chrom is not None

    if logger:
        logger.info('Loading sumstats')
    # Load data
    sumstats = fm_utils.load_sumstats(
        in_pq,
        study_id,
        phenotype_id=phenotype_id,
        biofeature=biofeature,
        chrom=chrom,
        min_maf=analysis_config['min_maf'],
        excl_mhc=analysis_config.get('exclude_MHC', None),
        logger=logger
    )
    if logger:
        logger.info('Completed loading: {0} rows'.format(sumstats.shape[0]))

    # Extract top loci
    if logger:
        logger.info('Starting top loci detection')
    top_loci = fm_top_loci.detect_top_loci(
        sumstats,
        in_plink,
        tmp_dir,
        method=method,
        maf=analysis_config['min_maf'],
        cojo_p=float(pval_threshold),
        cojo_window=analysis_config['cojo_wind'],
        cojo_collinear=analysis_config['cojo_colin'],
        clump_dist=analysis_config['clump_dist'],
        clump_p=float(pval_threshold),
        logger=logger
    )
    if logger:
        logger.info(
        'Completed top loci: {0} loci'.format(top_loci.shape[0]))

    # DEBUG only run for top loci
    # n_head = 3
    # if logger:
    #     logger.warning('Only running for top {0} loci'.format(n_head))
    # top_loci = top_loci.head(n_head)

    # Perform credible set analysis on each detected locus
    if logger:
        logger.info('Starting credible set analysis')
    credset_res_list = []
    for index_info in top_loci.to_dict(orient='records'):
        # Run
        credset_res = fm_credible_set.run_credible_set_for_locus(
            index_info,
            sumstats,
            top_loci,
            in_plink,
            tmp_dir,
            fm_wind=analysis_config['fm_wind'],
            cojo_window=analysis_config['cojo_wind'],
            cojo_collinear=analysis_config['cojo_colin'],
            pp_threshold=analysis_config['pp_threshold'],
            method=method,
            logger=logger)
        # Append result
        if credset_res is not None:
            credset_res_list.append(credset_res)

    # Concat credible sets together if they exist
    if len(credset_res_list) > 0:
        credset_results = pd.concat(credset_res_list)
    # Else creat an empty df
    else:
        credset_results = df_empty(
            columns=fm_utils.get_meta_info(type='cred_set').keys(),
            dtypes=fm_utils.get_meta_info(type='cred_set').values()
        )

    if logger:
        logger.info('Completed credible set analysis')

    return top_loci, credset_results

def df_empty(columns, dtypes, index=None):
    ''' Creat an empty df
    '''
    assert len(columns)==len(dtypes)
    df = pd.DataFrame(index=index)
    for c, d in zip(columns, dtypes):
        df[c] = pd.Series(dtype=d)
    return df
