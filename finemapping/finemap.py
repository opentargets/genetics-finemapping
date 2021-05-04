#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#
# Code for running FINEMAP (Benner et al.)

import utils as fm_utils
import subprocess as sp
import os
import pandas as pd
import traceback


def run_finemap_for_locus(
            index_info,
            sumstats,
            top_loci,
            in_plink,
            temp_dir,
            finemap_run_list,
            fm_wind,
            cojo_window,
            pp_threshold,
            logger=None):
    ''' Run FINEMAP at a given locus (specified by index_info)
    Args:
        index_info (dict): index locus information
        sumstats (pd.df): full summary stats
        top_loci (pd.df): table of top loci (for conditional analysis)
        finemap_run_list: 
    '''
    try:
        if logger:
            logger.info(
            '\n- Starting FINEMAP analysis at index variant {0}'.format(index_info['variant_id']))

        #temp_dir = os.path.join(temp_dir, 'finemap')
        chrom = str(sumstats.head(1)['chrom'].values[0])
        in_plink = in_plink.format(chrom=chrom)

        # Extract `cojo_window` region surrounding the index variant
        sumstat_wind = fm_utils.extract_window(
            sumstats, index_info['chrom'], index_info['pos'], cojo_window)
        if logger:
            logger.info('  {0} variants in {1}kb cojo window around index'
                        ' variant'.format(sumstat_wind.shape[0], cojo_window))

        # Get list of 'top loci' variants in the current window
        window_top_vars = tuple( set(sumstat_wind.variant_id).intersection(set(top_loci.variant_id)) )

        # Check if the "window top vars" are the same as a set we have previously done
        # if window_top_vars in finemap_res_list:
        #     logger.info('  Set of top variants {0} has already been done with FINEMAP. Skipping this set.'
        #                 .format(window_top_vars))
        #     return None
        if index_info['variant_id'] in finemap_run_list:
            logger.info('  Signal {0} has already been done with FINEMAP. Skipping this run.'
                        .format(index_info['variant_id']))
            return None
        
        # Extract `fm_wind` region surrounding the index variant
        sumstat_wind = fm_utils.extract_window(
            sumstats, index_info['chrom'], index_info['pos'], fm_wind)
        if logger:
            logger.info('  {0} variants in {1}kb fine-mapping window around index'
                        ' variant'.format(sumstat_wind.shape[0], fm_wind))

        #locus_name = ".".join(str(x).replace(':', '_') for x in sorted(window_top_vars))
        locus_name = index_info['variant_id']
        file_pref = make_file_name_prefix(sumstats.head(1))
        # Root name of files for finemap
        file_root = os.path.join(temp_dir, '{0}.{1}'.format(file_pref, locus_name))

        finemap_res = None
        if logger:
            logger.info('  Writing FINEMAP input files...')
        write_finemap_input(sumstat_wind, in_plink, file_root, pp_threshold=pp_threshold, logger=logger)

        if logger:
            logger.info('  Running FINEMAP for locus {0}'.format(locus_name))
        finemap_res = run_finemap(file_root, index_info['study_id'], locus_name, logger=logger)

        # Clean up large LD file immediately
        try:
            os.remove(file_root + '.ld')
        except OSError:
            pass

    except Exception as e:
        if logger:
            logger.info('  ERROR running FINEMAP...')
            logger.info(traceback.format_exc())
            logger.info(str(e))
        raise
    
    return locus_name, finemap_res

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


def run_finemap(file_root, study_id, locus_name, logger=None):
    # Constuct command
    finemap_corr_config = 0.9
    finemap_maxcausal_snps = 10
    cmd =  [
        'finemap --sss --log'.format(),
        '--in-files {0}'.format(file_root + '.master'),
        '--corr-config {0}'.format(finemap_corr_config),
        '--n-causal-snps {0}'.format(finemap_maxcausal_snps)
    ]

    # Run command
    # TODO: Capture output in a better way
    if logger:
        logger.info('  {0}'.format(' '.join(cmd)))
    fnull = open(os.devnull, 'w')
    cp = sp.run(' '.join(cmd), shell=True, stderr=sp.STDOUT)

    # Log error if return code is not 0
    if cp.returncode != 0:
        error = read_error_from_log(file_root + '.log')
        if logger:
            logger.error('  FINEMAP error:\n\n{0}\n'.format(error))
    
    finemap_snp_file = '{0}.snp'.format(file_root)
    if os.path.exists(finemap_snp_file):
        finemap_snps = pd.read_csv(finemap_snp_file, sep=' ', header=0)
        finemap_snps.insert(0, 'locus_name', [locus_name] * finemap_snps.shape[0], allow_duplicates = False)
        finemap_snps.insert(0, 'study_id', [study_id] * finemap_snps.shape[0], allow_duplicates = False)
    else:
        if logger:
            logger.warning(
                '  FINEMAP output not found. Returning empty df in place')
        finemap_snps = fm_utils.df_empty(
            columns=fm_utils.get_meta_info(type='finemap_snp').keys(),
            dtypes=fm_utils.get_meta_info(type='finemap_snp').values()
        )
    return(finemap_snps)


def write_finemap_input(sumstats, in_plink, file_root, pp_threshold, logger=None):
    ''' Writes files used by FINEMAP: <locus>.z, <locus>.ld
    '''
    finemap_master = file_root + '.master'
    finemap_z = file_root + '.z'
    finemap_ld = file_root + '.ld'
    # Output files
    finemap_snp = file_root + '.snp'
    finemap_config = file_root + '.config'
    finemap_cred = file_root + '.cred'
    finemap_log = file_root + '.log'

    #z;ld;snp;config;cred;log;n_samples
    #finemap/input/11_85867875.IGAP1_GWAX.z;finemap/input/11_85867875.IGAP1_GWAX.ld;finemap/out_ncausal_2/11_85867875.IGAP1_GWAX.2.snp;finemap/out_ncausal_2/11_85867875.IGAP1_GWAX.2.config;finemap/out_ncausal_2/11_85867875.IGAP1_GWAX.2.cred;finemap/out_ncausal_2/11_85867875.IGAP1_GWAX.2.log;166860
    with open(finemap_master, "w") as f:
        f.write("z;ld;snp;config;cred;log;n_samples\n")
        f.write(";".join([finemap_z, finemap_ld, finemap_snp, finemap_config, finemap_cred, finemap_log, str(10000)]) + "\n")
    
    # FINEMAP takes sumstats and LD as input. We use plink to compute LD, which
    # outputs the LD matrix in the same order as SNPs in the *.bim file - not in
    # the order specified in the snplist file passed in. So we must reorder the
    # sumstats to match the order in the LD file.

    # Read the genotype (.bim) file and subset to SNPs found there, and in the
    # same order
    bim = pd.read_csv(in_plink + ".bim", sep = '\t', header = None, names = ['chr', 'variant_id', 'cm', 'pos', 'a1', 'a2'])
    bim = bim.drop_duplicates(subset = 'variant_id')
    sumstats = sumstats.drop_duplicates(subset = 'variant_id')
    shared_variants = bim[['variant_id']].merge(sumstats['variant_id'], how="inner", on="variant_id")
    sumstats.index = sumstats['variant_id']
    sumstats_flt = sumstats.loc[shared_variants['variant_id'],:]
    if logger:
        logger.info('  {0} variants shared between LD reference and sumstats'.format(sumstats_flt.shape[0]))

    # Check that alleles match at all SNPs between genotype file and sumstats file
    bim.index = bim['variant_id']
    bim_flt = bim.loc[shared_variants['variant_id'],:]
    alleles_match = bim_flt.a1 == sumstats_flt.alt
    if logger:
        logger.info('  {0} of {1} variants have matching alleles'.format(sum(alleles_match), sumstats_flt.shape[0]))

    # Write sumstats to file for FINEMAP
    sumstat_to_finemap(sumstats_flt, finemap_z)

    # Write SNP list to use for plink to compute LD
    snplist_file = file_root + ".snplist"
    sumstats_flt.variant_id.to_csv(snplist_file, index=None, header=False)

    # Use plink to compute LD matrix
    # Constuct command
    cmd =  [
        'plink --bfile {0}'.format(in_plink),
        '--make-bed --keep-allele-order',
        '--extract {0}'.format(snplist_file),
        '--allow-extra-chr',
        '--out {0}'.format(file_root)
    ]
    # Command I ran in AD-finemap
    #os.system("plink --bfile gcta/input/ukbb_sample.{}.merged --extract {} --r square --out {}".format(chrom, snplist_fname, locus_fname))

    if logger:
        logger.info('  Running plink to extract SNPs: {0}'.format(' '.join(cmd)))

    # Run command
    fnull = open(os.devnull, 'w')
    plink_cp = sp.run(' '.join(cmd), shell=True, stdout=fnull, stderr=sp.STDOUT)

    cmd =  [
        'plink --bfile {0}'.format(in_plink),
        '--r square --keep-allele-order',
        '--extract {0}'.format(snplist_file),
        '--allow-extra-chr',
        '--out {0}'.format(file_root)
    ]
    if logger:
        logger.info('  Running plink to get LD: {0}'.format(' '.join(cmd)))

    # Run command
    fnull = open(os.devnull, 'w')
    plink_cp = sp.run(' '.join(cmd), shell=True, stdout=fnull, stderr=sp.STDOUT)

    # Log error if plink return code is not 0
    if plink_cp.returncode != 0:
        error = read_error_from_log(file_root + '.log')
        if logger:
            logger.error('  plink error:\n\n{0}\n'.format(error))
    
    # Convert tabs to spaces in LD output, necessary for FINEMAP
    cmd = "sed -i 's/\t/ /g' {0}.ld".format(file_root)
    if logger:
        logger.info('  Running sed to convert tabs to spaces'.format(cmd))
    sed_cp = sp.run(cmd, shell=True, stdout=fnull, stderr=sp.STDOUT)
    if sed_cp.returncode != 0 and logger:
        logger.error('  sed error converting tabs to spaces in LD file')

    
def sumstat_to_finemap(sumstats, outf):
    ''' Writes a sumstat df as a GCTA compliant file
    Args:
        sumstats (pd.df)
        outf (str): location to write to
        min_p (float): only write rows with p < p_threshold
    '''
    # Make temp dir if it doesn't exist
    os.makedirs(os.path.split(outf)[0], exist_ok=True)

    # Rename and extract required columns
    outdata = sumstats.rename(
        columns={"variant_id":"rsid",
                 "chrom":"chromosome",
                 "pos":"position",
                 "alt":"allele1",
                 "ref":"allele2",
                 "pval":"p"})
    outdata['maf'] = outdata['eaf'].apply(lambda f: f if f <= 0.5 else 1 - f)
    outdata = outdata.loc[:, ["rsid", "chromosome", "position", "allele1", "allele2", "maf", "beta", "se"]]

    # Save sumstats
    outdata.to_csv(outf, sep=' ', index=None)

    return 0


def read_error_from_log(log_file):
    ''' Reads a FINEMAP log file and returns lines containing the word "error"
    '''
    error_lines = []
    with open(log_file, 'r') as in_h:
        for line in in_h:
            if 'error' in line.lower():
                error_lines.append(line.rstrip())
    return '\n'.join(error_lines)


def make_file_name_prefix(row):
    ''' Takes a row of the sumstat dataframe and creates a unique file prefix
    Args:
        row (dask Series)
    Return:
        str
    '''
    cols = ['study_id', 'phenotype_id', 'bio_feature', 'chrom']
    pref = '_'.join([str(x) for x in row[cols].values[0]])
    return pref
