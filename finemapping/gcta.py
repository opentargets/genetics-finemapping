#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import subprocess as sp
import os
import pandas as pd

def get_conditional_top_loci(
        sumstats,
        in_plink,
        temp_dir,
        maf=0.01,
        cojo_p=5e-8,
        cojo_window=500,
        cojo_collinear=0.9,
        logger=None,
        split_ld=False
    ):
    ''' Uses GCTA-cojo to perform conditional analysis
    Args:

    '''

    # Extract filename prefix and chrom number
    chrom = str(sumstats.head(1)['chrom'].values[0])
    file_pref = make_file_name_prefix(sumstats.head(1))
    gcta_in = os.path.join(temp_dir, '{0}.gcta_format.tsv'.format(file_pref))
    gcta_snplist = os.path.join(temp_dir, '{0}.snplist.txt'.format(file_pref))
    gcta_out = os.path.join(temp_dir, '{0}.gcta_out'.format(file_pref))

    # Write sumstats
    sumstat_to_gcta(sumstats, gcta_in, gcta_snplist, p_threshold=cojo_p)

    ld_file = in_plink.format(chrom=chrom)
    if split_ld:
        # Determine if we can use a split LD reference, which speeds up GCTA dramatically.
        minpos = min(sumstats['pos'])
        maxpos = max(sumstats['pos'])
        minposMb = max(1, int(minpos / 1e6))
        maxposMb = max(1, int(maxpos / 1e6))
        if maxposMb - minposMb > 2:
            logger.error('split_ld parameter is true, but sumstat window is too large ({0} - {1}) to use split LD. Reverting to full chromosome LD.'.format(minpos, maxpos))
            split_ld = False
        else:
            # Take midpoint of window, in Mb
            MB_pos = int((maxposMb + minposMb) / 2)
            # We assume that the LD file per chromosome has been split into 3-Mb chunks
            window_size = int(3e6)
            ld_file = get_ld_fname(ld_file, MB_pos, window_size)

    # Construct command
    cmd = [
           'gcta64 --bfile {0}'.format(ld_file),
           '--chr {0}'.format(chrom),
           '--extract {0}'.format(gcta_snplist),
           '--maf {0}'.format(maf),
           '--cojo-p {0}'.format(cojo_p),
           '--cojo-wind {0}'.format(cojo_window),
           '--cojo-collinear {0}'.format(cojo_collinear),
           '--cojo-file {0}'.format(gcta_in),
           '--cojo-slct',
           '--out {0}'.format(gcta_out)
           ]

    if logger:
        logger.info('GCTA command:\n' + ' '.join(cmd))

    # Run command
    fnull = open(os.devnull, 'w')
    cp = sp.run(' '.join(cmd), shell=True, stdout=fnull, stderr=sp.STDOUT)

    # Log error if GCTA return code is not 0
    if cp.returncode != 0:
        gcta_log_file = '{0}.log'.format(gcta_out)
        gcta_error = read_error_from_gcta_log(gcta_log_file)
        logger.error('GCTA error:\n\n{0}\n'.format(gcta_error))

    # Read results (if they exist)
    gcta_res = '{0}.jma.cojo'.format(gcta_out)
    if os.path.exists(gcta_res):
        selected_snps = pd.read_csv(gcta_res, sep='\t', header=0)['SNP']
    else:
        selected_snps = []

    # Filter and return selected variants
    top_loci = sumstats.loc[sumstats['variant_id'].isin(selected_snps), :]

    return top_loci

def get_ld_fname(ld_file, MB_pos, window_size):
    ''' Determine the "split" LD file name, given the base (chromosome-level) name of an LD file,
        and check that it exists.
    '''
    def make_ld_fname(basefname, MB_pos):
        window_start = int(MB_pos * 1e6 - 1e6)
        window_end = int(window_start + window_size)
        return(basefname + '.{:d}_{:d}'.format(window_start, window_end))

    ld_fname = make_ld_fname(ld_file, MB_pos)
    # Check that the LD file exists
    if not os.path.exists(ld_fname + ".bim"):
        # If not, try the previous window, as we may be near the chromosome end
        ld_fname = make_ld_fname(ld_file, MB_pos - 1)
        if not os.path.exists(ld_fname + ".bim"):
            ld_fname = make_ld_fname(ld_file, MB_pos - 2)
    return ld_fname

def read_error_from_gcta_log(log_file):
    ''' Reads a GCTA log file and returns lines containing the word "error"
    '''
    error_lines = []
    with open(log_file, 'r') as in_h:
        for line in in_h:
            if 'error' in line.lower():
                error_lines.append(line.rstrip())
    return '\n'.join(error_lines)

def perform_conditional_adjustment(sumstats,
        in_plink,
        temp_dir,
        index_var,
        chrom,
        condition_on,
        cojo_window=500,
        cojo_collinear=0.9,
        logger=None,
        split_ld=False,
        var_pos=None):
    ''' Uses GCTA-cojo to perform conditional analysis
    Args:
        sumstats (pd.df)
        in_plink (str)
        temp_dir (str)
        index_var (str)
        chrom (str)
        condition_on (list): list of variants to condition on
        split_ld (bool): whether to use split LD files or not
        var_pos (int): needed if split_ld is True: integer position of index_var
    Return:
        pd.df of conditionally adjusted summary stats
    '''
    # Extract filename prefix and chrom number
    file_pref = make_file_name_prefix(sumstats.head(1))
    gcta_in = os.path.join(temp_dir, '{0}.{1}.gcta_format.tsv'.format(
        file_pref,
        index_var.replace(':', '_')))
    gcta_snplist = os.path.join(temp_dir, '{0}.{1}.snplist.txt'.format(
        file_pref,
        index_var.replace(':', '_')))
    gcta_out = os.path.join(temp_dir, '{0}.{1}.gcta_out'.format(
        file_pref,
        index_var.replace(':', '_')))
    gcta_cond = os.path.join(temp_dir, '{0}.{1}.cond_list.txt'.format(
        file_pref,
        index_var.replace(':', '_')))

    # Write sumstats
    sumstat_to_gcta(sumstats, gcta_in, gcta_snplist)

    # Write a conditional list
    write_cond_list(condition_on, gcta_cond)

    ld_file = in_plink.format(chrom=chrom)
    if split_ld:
        # We assume that the LD file per chromosome has been split into 3-Mb chunks
        window_size = int(3e6)
        MB_pos = max(1, int(var_pos / 1e6))
        ld_file = get_ld_fname(ld_file, MB_pos, window_size)

    # Constuct command
    cmd =  [
        'gcta64 --bfile {0}'.format(ld_file),
        '--chr {0}'.format(chrom),
        '--extract {0}'.format(gcta_snplist),
        '--cojo-file {0}'.format(gcta_in),
        '--cojo-cond {0}'.format(gcta_cond),
        '--out {0}'.format(gcta_out)
    ]

    if logger:
        logger.info('GCTA command:\n' + ' '.join(cmd))

    # Run command
    fnull = open(os.devnull, 'w')
    cp = sp.run(' '.join(cmd), shell=True, stdout=fnull, stderr=sp.STDOUT)

    # Log error if GCTA return code is not 0
    if cp.returncode != 0:
        gcta_log_file = '{0}.log'.format(gcta_out)
        gcta_error = read_error_from_gcta_log(gcta_log_file)
        if logger:
            logger.error('  GCTA error:\n\n{0}\n'.format(gcta_error))

    # Read output if it exists
    gcta_res = '{0}.cma.cojo'.format(gcta_out)
    if os.path.exists(gcta_res):
        cond_res = pd.read_csv(gcta_res, sep='\t', header=0)
        # Merge to sumstats
        sumstat_cond = merge_conditional_w_sumstats(sumstats, cond_res)
    # Otherwise return empty sumstat df with conditional columns added
    else:
        if logger:
            logger.warning(
                '  GCTA output not found. Returning empty df in place'.format(gcta_res))
        sumstat_cond = sumstats.loc[[], :]
        for col in ['beta_cond', 'se_cond', 'pval_cond']:
            sumstat_cond[col] = None

    return sumstat_cond

def merge_conditional_w_sumstats(sumstats, cond_res):
    ''' Adds the results of the conditional analysis as columns to the
        summary stats
    Args:
        sumstats (pd.df)
        cond_res (pd.df)
    Returns:
        pd.df of sumstats
    '''
    # Make a variant_id key on the conditional results
    cond_res = cond_res.rename(columns={'SNP': 'variant_id'})
    
    # Merge
    merged = pd.merge(
        sumstats,
        cond_res.loc[:, ['variant_id', 'bC', 'bC_se', 'pC']],
        on='variant_id',
        how='inner'
    )
    return merged.rename(columns={
        'bC':'beta_cond',
        'bC_se':'se_cond',
        'pC':'pval_cond',
    })

def write_cond_list(cond_list, outf):
    ''' Write a list of variants to a file
    '''
    with open(outf, 'w') as out_h:
        for var in cond_list:
            out_h.write(var.replace('_', ':') + '\n')
    return outf

def sumstat_to_gcta(sumstats, outf, snplist, p_threshold=None):
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
        columns={"variant_id":"SNP",
                 "alt":"A1",
                 "ref":"A2",
                 "eaf":"freq",
                 "beta":"b",
                 "pval":"p",
                 "n_total":"N"})
    outdata = outdata.loc[:, ["SNP", "A1", "A2", "freq", "b", "se", "p", "N"]]

    # Remove rows where p > p_threshold
    if p_threshold:
        outdata = outdata.loc[outdata['p'] <= p_threshold, :]

    # Save sumstats
    outdata.to_csv(outf, sep='\t', index=None)

    # Save snplist
    outdata.SNP.to_csv(snplist, index=None, header=False)

    return 0

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
