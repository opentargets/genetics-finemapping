#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import finemapping.utils
import finemapping.gcta
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
            method='conditional'):
    ''' Run credible set analysis at a given locus (speficied by index_info)
    Args:
        index_info (dict): index locus information
        sumstats (pd.df): full summary stats
        top_loci (pd.df): table of top loci (for conditional analysis)
        method ([conditional|distance]): whether to perform conditional analysis
            or extract a set distance
    '''

    temp_dir = os.path.join(temp_dir, 'credible_set')

    # Extract region surrounding the index variant
    sumstat_wind = finemapping.utils.extract_window(
        sumstats, index_info['chrom'], index_info['pos'], fm_wind
    )

    # Perform conditional analysis
    sumstat_cond = None
    if method == 'conditional':

        # Get list of variants to condition on
        cond_list = make_list_to_condition_on(
            index_info['variant_id'],
            sumstat_wind.variant_id,
            top_loci.variant_id
        )

        # Only do conditional if there are variants to condition on
        if len(cond_list) > 0:
            sumstat_cond = finemapping.gcta.perfrom_conditional_adjustment(
                sumstat_wind,
                in_plink,
                temp_dir,
                index_info['variant_id'],
                index_info['chrom'],
                cond_list
            )

    # If conditional analysis was not perform, we need to make the sumstat df
    # look the same, add (beta_cond, se_cond, pval_cond) columns
    if sumstat_cond is None:
        sumstat_cond = sumstat_wind
        sumstat_cond['beta_cond'] = sumstat_cond['beta']
        sumstat_cond['se_cond'] = sumstat_cond['se']
        sumstat_cond['pval_cond'] = sumstat_cond['pval']

    # Do credible set analysis
    cred_sets = calc_credible_sets(sumstat_cond)

    # Add index variant as a column
    cred_sets.loc[:, 'index_variant_id'] = index_info['variant_id']

    return cred_sets

def calc_credible_sets(data):
    ''' Calculates credible sets from provided sumstats
    Args:
        data (pd.df): sumstats to perform analysis on
    Returns
        pd.df of results
    '''
    # Need to fix pC == 0.0 as this will produce inf logABF. Set to sys.float_info.min
    if (data.pval_cond == 0.0).any():
        print("Warning: some pval_cond == 0.0 in {0}\n - setting to sys.float_info.min")
        data.pval_cond[data.pval_cond == 0.0] = sys.float_info.min

    # Calculate case proportions
    data['case_prop'] = data['n_cases'] / data['n_samples']

    # Calc ABFs
    # print(calc_abf(0.808621, 0.17690, 290365, 0.6203537)) # Should return -3.311501
    data["logABF"] = data.apply(
        lambda row: calc_abf(pval=row['pval_cond'],
                             maf=freq_to_maf(row['eaf']),
                             n=row['n_samples'],
                             prop_cases=row['case_prop'] if row['is_cc'] else None
        ), axis=1)
    data = data.sort_values("logABF", ascending=False)

    # Calculate posterior probability for each SNP
    sum_lABF = log_sum(data["logABF"])
    data["postprob"] = data["logABF"].apply(np.exp) / np.exp(sum_lABF)

    # Calc cumulative sum of the posterior probabilities
    data["postprob_cumsum"] = data["postprob"].cumsum()

    # Find 99% and 95% credible sets - this is horrible
    set_idx = data["postprob_cumsum"].gt(0.99).tolist().index(True)
    data["is99_credset"] = [1] * (set_idx + 1) + [0] * (data.shape[0] - (set_idx + 1))
    set_idx = data["postprob_cumsum"].gt(0.95).tolist().index(True)
    data["is95_credset"] = [1] * (set_idx + 1) + [0] * (data.shape[0] - (set_idx + 1))

    # Only keep rows that are in the 95 or 99% credible sets
    to_keep = ((data["is95_credset"] == 1) | (data["is99_credset"] == 1))
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
