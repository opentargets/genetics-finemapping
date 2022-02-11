#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import utils as fm_utils
import gcta as fm_gcta
import os
import numpy as np
import pandas as pd
from scipy.stats import norm
import sys

def run_credible_set_for_locus(
            index_info,
            sumstats,
            top_loci,
            in_plink,
            temp_dir,
            fm_wind,
            cojo_window,
            cojo_collinear,
            pp_threshold,
            method='conditional',
            logger=None,
            split_ld=False):
    ''' Run credible set analysis at a given locus (speficied by index_info)
    Args:
        index_info (dict): index locus information
        sumstats (pd.df): full summary stats
        top_loci (pd.df): table of top loci (for conditional analysis)
        method ([conditional|distance|none]): whether to perform conditional analysis
            or extract a set distance
    '''

    if logger:
        logger.info(
        '\n- Running credible set analysis for {0}'.format(index_info['variant_id']))

    temp_dir = os.path.join(temp_dir, 'credible_set')

    # Perform conditional analysis
    sumstat_cond = None
    if method == 'conditional':
        # Extract `cojo_window` region surrounding the index variant
        sumstat_wind = fm_utils.extract_window(
            sumstats, index_info['chrom'], index_info['pos'], cojo_window)
        if logger:
            logger.info('  {0} variants in {1}kb cojo window around index'
                        ' variant'.format(sumstat_wind.shape[0], cojo_window))

        # Get list of variants to condition on
        cond_list = make_list_to_condition_on(
            index_info['variant_id'],
            sumstat_wind.variant_id,
            top_loci.variant_id
        )

        if logger:
            logger.info('  conditioning on {0} variants'.format(
            len(cond_list)))

        # Only do conditional if there are variants to condition on
        if len(cond_list) > 0:
            sumstat_cond = fm_gcta.perform_conditional_adjustment(
                sumstat_wind,
                in_plink,
                temp_dir,
                index_info['variant_id'],
                index_info['chrom'],
                cond_list,
                cojo_window,
                cojo_collinear,
                logger=logger,
                split_ld=split_ld,
                var_pos=index_info['pos']
            )
    else:
        sumstat_wind = sumstats
        if logger:
            logger.info('  not conditioning as method != conditional')

    # If conditional analysis was not performed, we need to make the sumstat df
    # look the same, add (beta_cond, se_cond, pval_cond) columns
    if sumstat_cond is None:
        sumstat_cond = sumstat_wind
        sumstat_cond['beta_cond'] = sumstat_cond['beta']
        sumstat_cond['se_cond'] = sumstat_cond['se']
        sumstat_cond['pval_cond'] = sumstat_cond['pval']

    # Extract `fm_wind` region surrounding the index variant
    # TODO: test results compared to Ed's code before, which didn't use sumstat_cond here
    sumstat_wind = fm_utils.extract_window(
        sumstat_cond, index_info['chrom'], index_info['pos'], fm_wind)
    if logger:
        logger.info('  {0} variants in {1}kb fine-mapping window around index'
                    ' variant'.format(sumstat_wind.shape[0], fm_wind))

    # Do credible set analysis
    if sumstat_wind.shape[0] > 0:

        if logger:
            logger.info('  calculating credible sets...')
        cred_sets = calc_credible_sets(sumstat_wind, pp_threshold=pp_threshold)

        if logger:
            logger.info('  found {0} in 95% and {1} in 99% cred sets'.format(
            cred_sets.is95_credset.sum(), cred_sets.is99_credset.sum()
            ))
            logger.info('  kept {0} vars with PP > {1}'.format(
            cred_sets.shape[0], pp_threshold
            ))

        # Script will fail if cred_sets is empty
        if cred_sets.shape[0] > 0:

            # Add index variant columns
            cred_sets.loc[:, 'lead_variant_id'] = index_info['variant_id']
            cred_sets[['lead_chrom', 'lead_pos', 'lead_ref', 'lead_alt']] = \
                cred_sets.lead_variant_id.str.split(':', expand=True)

            # Add column specifying method used
            cred_sets.loc[:, 'multisignal_method'] = method

            # Format output table
            cred_sets = format_credset_output(cred_sets)
        
        # Else None if credible set results is empty
        else:
            cred_sets = None

    # If df is empty skip analysis
    else:
        if logger:
            logger.warning('  skipping credible set analysis')
        cred_sets = None

    return cred_sets

def format_credset_output(cred_sets):
    ''' Formats the cred_sets table for output
    Args:
        cred_sets (pd.df)
    Returns:
        pd.df
    '''
    cols = fm_utils.get_credset_out_columns()
    meta = fm_utils.get_meta_info(type='cred_set')
    df = (
        cred_sets.loc[:, cols.keys()]
                 .rename(columns=cols)
                 .astype(dtype=meta)
    )
    return df

def calc_credible_sets(data, pp_threshold):
    ''' Calculates credible sets from provided sumstats
    Args:
        data (pd.df): sumstats to perform analysis on
        pp_threshold (float): returns any variant in ( (95% OR 99% threshold) OR pp > pp_threshold )
    Returns
        pd.df of results
    '''
    # Need to fix pC == 0.0 as this will produce inf logABF. Set to sys.float_info.min
    if (data.pval_cond == 0.0).any():
        print("Warning: some pval_cond == 0.0 in {0}\n - setting to sys.float_info.min")
        data.pval_cond[data.pval_cond == 0.0] = sys.float_info.min

    # Calculate case proportions
    data['case_prop'] = data['n_cases'] / data['n_total']

    # Calc ABFs
    # print(calc_abf(0.808621, 0.17690, 290365, 0.6203537)) # Should return -3.311501
    data["logABF"] = data.apply(
        lambda row: calc_abf(pval=row['pval_cond'],
                             maf=freq_to_maf(row['eaf']),
                             n=row['n_total'],
                             prop_cases=row['case_prop'] if row['is_cc'] else None
        ), axis=1)
    data = data.sort_values("logABF", ascending=False)

    # Calculate posterior probability for each SNP
    sum_lABF = log_sum(data["logABF"])
    data["postprob"] = (data["logABF"] - sum_lABF).apply(np.exp)

    # Calc cumulative sum of the posterior probabilities
    data["postprob_cumsum"] = data["postprob"].cumsum()

    # Find 99% and 95% credible sets - this is horrible
    set_idx = data["postprob_cumsum"].gt(0.95).tolist().index(True)
    data["is95_credset"] = [1] * (set_idx + 1) + [0] * (data.shape[0] - (set_idx + 1))
    data["is95_credset"] = data["is95_credset"].map({1:True, 0:False})
    set_idx = data["postprob_cumsum"].gt(0.99).tolist().index(True)
    data["is99_credset"] = [1] * (set_idx + 1) + [0] * (data.shape[0] - (set_idx + 1))
    data["is99_credset"] = data["is99_credset"].map({1:True, 0:False})

    # Only keep rows that are in the 95 or 99% credible sets
    to_keep = ((data["is95_credset"] | data["is99_credset"])
               & (data["postprob"] > pp_threshold) )
    cred_set_res = data.loc[to_keep, :]

    return cred_set_res

def make_list_to_condition_on(index_var, all_vars, top_vars):
    ''' Makes a list of variants on which to condition on
    Args:
        index_var (str): index variant at this locus
        all_vars (pd.Series): A series of all variant IDs in the locus window
        top_vars (pd.Series): A series of top loci variant IDs
    Returns:
        list of variants to condition on
    '''
    window_top_vars = set(all_vars).intersection(set(top_vars))
    cond_list = list(window_top_vars - set([index_var]))
    return cond_list

def calc_abf(pval, maf, n, prop_cases=None):
    """ Caluclate Approximate Bayes Factor (Wakefield, 2009, Genet Epidemiol.).
        Based on code from coloc: https://github.com/chr1swallace/coloc
    Args:
        pval (float): GWAS p-value
        maf (float): Minor allele freq
        n (int): Sample size
        prop_cases (float or None): number of cases, if left blank will assume
            quantitative trait
    Returns:
        natural log(ABF)
    """
    # Assert/set types
    pval = float(pval)
    maf = float(maf)
    n = int(n)
    prop_cases = float(prop_cases) if prop_cases else None

    # Estimate variance for quant trait
    if prop_cases is None:
        sd_prior = 0.15
        v = var_data(maf, n)
    # Estimate var for cc study
    else:
        sd_prior = 0.2
        v = var_data_cc(maf, n, prop_cases)

    # Calculate Z-score
    z = np.absolute(norm.ppf(pval / 2))

    # Calc shrinkage factor: ratio of the prior variance to the total variance
    r = sd_prior**2 / (sd_prior**2 + v)

    # Approximate BF - ln scale to compare in log natural scale with LR diff
    lABF = 0.5 * (np.log(1 - r) + (r * z**2))

    return lABF

def log_sum(l):
    """ Calculates the log of the sum of the exponentiated logs taking out the
        max, i.e. insuring that the sum is not Inf
    Args:
        l (pandas Series)
    Returns:
        Sum of exponentiated logs
    """
    l_max = l.max()
    l_logsum = l_max + np.log(np.sum(np.exp(l - l_max)))
    return l_logsum

def freq_to_maf(freq):
    """ Convert allele frequency to minor allele freq
    """
    return min(freq, 1-freq)

def var_data(maf, n):
    """ Calc variance of MLE of beta for quantitative trait, assuming var(y)=0
    """
    var = 1 / (2 * n * maf * (1 - maf))
    return var

def var_data_cc(maf, n, prop_cases):
    """ Calc variance of MLE of beta for case-control
    """
    var = 1 / (2 * n * maf * (1 - maf) * prop_cases * (1 - prop_cases))
    return var
