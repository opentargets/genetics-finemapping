#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Jeremy Schwartzentruber
#

import sys
import os
import argparse
import logging
import pprint
from datetime import datetime
import yaml
import pandas as pd

# Import utils file from finemapping folder
# https://stackoverflow.com/questions/4383571/importing-files-from-different-folder
sys.path.insert(1, '../finemapping')
import utils as fm_utils

def main():

    # Args
    args = parse_args()
    start_time = datetime.now()

    # Make output and temp directories
    os.makedirs(os.path.split(args.toploci)[0], exist_ok=True)
    os.makedirs(os.path.split(args.credset)[0], exist_ok=True)
    os.makedirs(os.path.split(args.log)[0], exist_ok=True)

    # Start logging
    logger = make_logger(args.log)
    logger.info('Started FinnGen finemapping ingest')
    logger.info('Printing command: \n' + ' '.join(sys.argv))
    logger.info('Printing args: \n' + pprint.pformat(vars(args), indent=2))

    # Load analysis config file
    with open(args.config_file, 'r') as in_h:
        config_dict = yaml.load(in_h)
    logger.info('Analysis config: \n' + pprint.pformat(config_dict, indent=2))
    
    # Read all SNPs into dataframe, then filter on SNP probability
    # (If this is inefficient, we could read the file in chunks and filter each chunk:
    # https://stackoverflow.com/questions/13651117/how-can-i-filter-lines-on-load-in-pandas-read-csv-function)
    snps = pd.read_csv(args.in_snp, sep="\t", compression="gzip")
    snps = snps.rename(columns={'prob': 'postprob'})
    if logger:
        logger.info('  {0} variants in fine-mapping window: {1}'.format(
            snps.shape[0], snps['region'][0]))

    # First remove SNPs without a credible set ID (they are labeled -1)
    # or with probability below our threshold
    snps = snps[snps['cs'] > 0]
    snps = snps[snps['postprob'] >= config_dict['pp_threshold']]
    if logger:
        logger.info('  kept {0} vars in credible sets and with PP > {1}'.format(
        snps.shape[0], config_dict['pp_threshold']
        ))
    if snps.shape[0] == 0:
        # There are sometimes no SNPs in credible sets. I don't know why this
        # happens in FinnGen's SuSIE output.
        if logger:
            logger.info('No SNPs in credible set - exiting.')
        return 0

    # Sort by SNP probability, and keep top SNP per locus
    snps_by_signal = snps.groupby(['region', 'cs'], as_index=False)

    def get_top_snp(df):
        top_snp = df.sort_values(by=['postprob'], ascending=False).head(1)
        if top_snp.loc[top_snp.index[0], 'p'] > float(config_dict['gwas_pval_threshold']):
            top_snp = top_snp.head(0)
        return(top_snp)
    top_loci_by_signal = snps_by_signal.apply(get_top_snp)

    if top_loci_by_signal.shape[0] == 0:
        # There are sometimes no SNPs with p < threshold in credible sets.
        if logger:
            logger.info('No SNPs in any credible set with raw p < {0} - exiting.'.format(float(config_dict['gwas_pval_threshold'])))
        return 0

    # Finally, since we are keeping only one SNP per region, keep only
    # the top SNP by p value per region.
    top_loci_by_region = top_loci_by_signal.groupby(['region'], as_index=False)
    def get_top_snp_by_pval(df):
        top_snp = df.sort_values(by=['p'], ascending=True).head(1)
        return(top_snp)
    top_loci = top_loci_by_region.apply(get_top_snp_by_pval)
    
    # Add additional required columns
    top_loci.loc[:, 'study_id'] = args.study_id
    top_loci.loc[:, 'clump_method'] = "finngen"
    top_loci.loc[:, 'phenotype_id'] = None
    top_loci.loc[:, 'bio_feature'] = None
    top_loci = top_loci.rename(columns={
        'v': 'variant_id',
        'chromosome': 'chrom',
        'position': 'pos',
        'allele1': 'ref',
        'allele2': 'alt',
        'p': 'pval'})
    top_loci['chrom'] = top_loci['chrom'].str.replace('chr', '')
    top_loci = format_top_loci_output(top_loci)
    # Needs to be added after format_... call (kind of a bug/hack)
    top_loci.insert(0, 'type', "gwas")

    # Compute 95% and 99% credset SNPs
    if logger:
        logger.info('  formatting credible sets...')

    cred_sets = snps_by_signal.apply(get_credible_sets, gwas_pval_threshold=float(config_dict['gwas_pval_threshold']))
    #cred_sets = []
    #for locus, frame in snps_by_signal:
    #    cred_set = frame.sort_values(by=['postprob'], ascending=False).apply(calc_credible_sets)
    #    cred_sets.append([cred_set])
    #red_sets = pd.DataFrame(cred_sets)

    # if logger:
    #     logger.info('  found {0} in 95% and {1} in 99% cred sets'.format(
    #     cred_sets.is95_credset.sum(), cred_sets.is99_credset.sum()
    #     ))

    # Add index variant columns
    cred_sets = cred_sets.groupby(['region', 'cs'], as_index=False).apply(add_lead_variant_cols)

    # Add additional required columns
    cred_sets.loc[:, 'multisignal_method'] = "SuSIE"
    cred_sets.loc[:, 'study_id'] = args.study_id
    cred_sets.loc[:, 'phenotype_id'] = None
    cred_sets.loc[:, 'bio_feature'] = None

    # Rename columns to match credset results format
    cred_sets = cred_sets.rename(columns={
        'v': 'variant_id',
        'chromosome': 'chrom',
        'position': 'pos',
        'allele1': 'ref',
        'allele2': 'alt',
        'p': 'pval'})
    cred_sets['chrom'] = cred_sets['chrom'].str.replace('chr', '')

    # FinnGen didn't do conditional analysis, so just use the original beta, etc
    # I'm not sure if these fields matter.
    cred_sets.loc[:, 'beta_cond'] = cred_sets.loc[:, 'beta']
    cred_sets.loc[:, 'se_cond'] = cred_sets.loc[:, 'se']
    cred_sets.loc[:, 'pval_cond'] = cred_sets.loc[:, 'pval']

    # I'm not sure if we need to calculate the logABF for FinnGen. It's part of the
    # credset table, but I wouldn't consider it a normal fine-mapping output. It's
    # only used in the conditional analysis / single-causal variant approach used
    # for our other GWAS.
    cred_sets['logABF'] = None

    # Format output table
    cred_sets = format_credset_output(cred_sets)

    # Needs to be added after format_credset_output call (kind of a bug/hack)
    cred_sets.insert(0, 'type', "gwas")

    # Write output
    logger.info('Writing outputs')

    # As json
    top_loci.to_json(
        args.toploci,
        orient='records',
        lines=True,
        compression='gzip',
        double_precision=15
    )
    cred_sets.to_json(
        args.credset,
        orient='records',
        lines=True,
        compression='gzip',
        double_precision=15
    )

    # # As parquet
    # top_loci.to_parquet(
    #     fname=args.toploci,
    #     engine='pyarrow',
    #     flavor='spark',
    #     compression='snappy',
    #     index=False
    # )
    # cred_sets.to_parquet(
    #     fname=args.credset,
    #     engine='pyarrow',
    #     flavor='spark',
    #     compression='snappy',
    #     index=False
    # )

    # Log time taken
    logger.info('Time taken: {0}'.format(
        datetime.now() - start_time
    ))
    logger.info('Finished!')

    return 0

def freq_to_maf(freq):
    """ Convert allele frequency to minor allele freq
    """
    return min(freq, 1-freq)

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

def format_top_loci_output(top_loci):
    ''' Formats the top_loci table for output
    Args:
        top_loci (pd.df)
    Returns:
        pd.df
    '''
    cols = fm_utils.get_toploci_out_columns()
    meta = fm_utils.get_meta_info(type='top_loci')
    df = (
        top_loci.loc[:, cols.keys()]
                .rename(columns=cols)
                .astype(dtype=meta)
    )
    return df

def get_credible_sets(data, gwas_pval_threshold):
    ''' Adds columns for credible set output. Assumes that the input data
        are already filtered for variants in 95% credset (as in FinnGen).
    Args:
        data (pd.df): dataframe with 'postprob' column
    Returns
        pd.df of results (only variants in 99% credible set)
    '''
    # Calc cumulative sum of the posterior probabilities
    data = data.sort_values(by=['postprob'], ascending=False)
    if data.loc[data.index[0], 'p'] > gwas_pval_threshold:
        data = data.head(0)
    else:
        data["postprob_cumsum"] = data["postprob"].cumsum()

        # Find 99% and 95% credible sets - this is horrible
        data["is95_credset"] = True
        data["is99_credset"] = True

        # Only keep rows that are in the 95 or 99% credible sets
        #to_keep = ((data["is95_credset"] | data["is99_credset"]))
        #cred_set_res = data.loc[to_keep, :]

    return data

def add_lead_variant_cols(df):
    # Data should already be sorted, but just in case...
    df = df.sort_values(by=['postprob'], ascending=False)
    leadvar = df.iloc[0,:]
    df.loc[:, 'lead_variant_id'] = leadvar['v']
    df[['lead_chrom', 'lead_pos', 'lead_ref', 'lead_alt']] = \
            leadvar['v'].split(':')
    cols = ['region', 'v', 'chromosome', 'position', 'allele1', 'allele2', 'maf', 'postprob', 'cs', 'postprob_cumsum', 'is95_credset', 'is99_credset', 'lead_variant_id', 'lead_chrom', 'lead_pos', 'lead_ref', 'lead_alt']
    return df

def parse_args():
    ''' Load command line args
    '''
    p = argparse.ArgumentParser()

    # Input files/params
    p.add_argument('--in_snp',
                   metavar="<file>",
                   help=("Input: SUSIE SNP file fine-mapping output"),
                   type=str,
                   required=True)
    p.add_argument('--config_file',
                   metavar="<file>",
                   help=("Input: analysis config file"),
                   type=str,
                   required=True)
    p.add_argument('--study_id',
                   metavar="<file>",
                   help=("Input: analysis config file"),
                   type=str,
                   required=True)
    # Output files
    p.add_argument('--toploci',
                   metavar="<file>",
                   help=("Output: top loci json file"),
                   type=str,
                   required=True)
    p.add_argument('--credset',
                   metavar="<file>",
                   help=("Output: credible set json file"),
                   type=str,
                   required=True)
    p.add_argument('--log',
                   metavar="<file>",
                   help=("Output: log file"),
                   type=str,
                   required=True)

    args = p.parse_args()

    # Convert "None" strings to None type
    for arg in vars(args):
        if getattr(args, arg) == "None":
            setattr(args, arg, None)

    return args

def make_logger(log_file):
    ''' Creates a logging handle.
    '''
    # Basic setup
    logging.basicConfig(
        level=logging.INFO,
        datefmt='%Y-%m-%d %H:%M:%S',
        stream=None)
    # Create formatter
    logFormatter = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
    rootLogger = logging.getLogger(__name__)
    # Add file logging
    fileHandler = logging.FileHandler(log_file, mode='w')
    fileHandler.setFormatter(logFormatter)
    rootLogger.addHandler(fileHandler)
    # Add stdout logging
    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    rootLogger.addHandler(consoleHandler)
     # Prevent logging from propagating to the root logger
    rootLogger.propagate = 0

    return rootLogger

if __name__ == '__main__':

    main()
