#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import finemapping.utils
import finemapping.top_loci
import finemapping.credible_set
import tempfile
import yaml
import json
import pandas as pd
from pprint import pprint # DEBUG
import dask
from dask.distributed import Client
from dask.distributed import Client
import os

def run_all_studies(gwas_sumstats_dir, mol_sumstats_dir, manifest_file,
        results_dir, analysis_config_file, temp_dir, skip_completed,
        skipped_results, n_cores=1):
    ''' Run the fine-mapping pipeline on all studies in the manifest file.
        Optionally, skips studies that already have results.
    Args:
        gwas_sumstats_dir (dir): directory containing gwas summary stats
        mol_sumstats_dir (dir): dir containing molecular trait gwas sumstats
        manifest_file (file): manifest file specifying which studies to run
        plink_ld (str): pattern matching ld reference for european samples
        results_dir (dir): where to save the results
        analysis_config_file (file): config file containing analysis parameters
        temp_dir (dir): where to store temporary files
        skip_completed (bool): whether studies with existing results should be
                               skipped
        skipped_results (dir): location of results for skipped studies
        n_cores (int): number of cores to use
    '''

    # Set scheduler params
    client = Client(n_workers=n_cores)

    # Load the manifest
    with open(manifest_file, 'r') as in_h:
        manifest = json.load(in_h)
    # DEBUG
    print('Manifest: ', end='')
    pprint(manifest)

    # Create delayed function
    run_single_study_delayed = dask.delayed(run_single_study, nout=2)

    # Initiate delayed results lists
    top_loci_results_list = []
    cred_set_results_list = []

    # Created delayed execution for each study in the manfiest
    for study_info in manifest:

        # Create pq input variable
        if study_info['type'] == 'gwas':
            in_pq = os.path.join(gwas_sumstats_dir, study_info['study_id'])
        elif study_info['type'] == 'molecular_qtl':
            in_pq = os.path.join(mol_sumstats_dir, study_info['study_id'])

        # Create temp input variable
        tmp_dir = os.path.join(temp_dir, study_info['study_id'])

        # Submit delayed
        top_loci, credset = run_single_study_delayed(
            in_pq=in_pq,
            in_plink=study_info['ld_ref'],
            study_id=study_info['study_id'],
            cell_id=study_info['cell_id'],
            gene_id=study_info['gene_id'],
            group_id=study_info['group_id'],
            chrom=study_info['chrom'],
            analysis_config_file=analysis_config_file,
            tmp_dir=tmp_dir,
            method=study_info['method']
        )

        # Add results to lists
        top_loci_results_list.append(top_loci)
        cred_set_results_list.append(credset)

    # Make dask.dfs from the delayed objects
    top_loci_dd = dask.dataframe.from_delayed(
        top_loci_results_list,
        meta=finemapping.utils.get_meta_info(type='top_loci'))
    cred_set_dd = dask.dataframe.from_delayed(
        cred_set_results_list,
        meta=finemapping.utils.get_meta_info(type='cred_set'))


    dask.compute(top_loci_dd, cred_set_dd)
    # top_loci_dd.compute()
    # cred_set_dd.compute()

    return 0

def run_single_study(in_pq,
                     in_plink,
                     study_id,
                     cell_id=None,
                     gene_id=None,
                     group_id=None,
                     chrom=None,
                     analysis_config_file=None,
                     tmp_dir=tempfile.gettempdir(),
                     method='conditional'):
    ''' Runs the top loci and credible set analysis on a single study
    '''

    pd.options.mode.chained_assignment = None # Silence pandas warning

    # TEMP it is reqired that chrom is not None
    assert chrom is not None

    # Load analysis config file
    with open(analysis_config_file, 'r') as in_h:
        analysis_config = yaml.load(in_h)
        # print('Config: ', analysis_config) # DEBUG

    # Load data
    sumstats = finemapping.utils.load_sumstats(
        in_pq,
        study_id,
        cell_id=cell_id,
        gene_id=gene_id,
        group_id=group_id,
        chrom=chrom,
        min_maf=analysis_config['min_maf'],
        excl_mhc=analysis_config.get('exclude_MHC', None),
        build=analysis_config['assembly']
    )

    # Extract top loci
    top_loci = finemapping.top_loci.detect_top_loci(
        sumstats,
        in_plink,
        tmp_dir,
        method=method,
        maf=analysis_config['min_maf'],
        cojo_p=float(analysis_config['cojo_p']),
        cojo_window=analysis_config['cojo_wind'],
        cojo_collinear=analysis_config['cojo_colin'],
        clump_dist=analysis_config['clump_dist'],
        clump_p=float(analysis_config['clump_pval'])
    )

    # DEBUG only run for top loci
    print('Warning: Only running for subset of loci')
    top_loci = top_loci.head(3)

    # Perform credible set analysis on each detected locus
    credset_res_list = []
    for index_info in top_loci.to_dict(orient='records'):
        # Run
        credset_res = finemapping.credible_set.run_credible_set_for_locus(
            index_info,
            sumstats,
            top_loci,
            in_plink,
            tmp_dir,
            fm_wind=analysis_config['fm_wind'],
            method=method
        )
        # Append result
        credset_res_list.append(credset_res)

    # Concat credible sets together if they exist
    if len(credset_res_list) > 0:
        credset_results = pd.concat(credset_res_list)
    # Otherwise create an empty df with matching meta data
    else:
        cols = finemapping.utils.get_credset_out_columns().values()
        meta = finemapping.utils.get_meta_info(type='cred_set')
        credset_results = pd.DataFrame(columns=cols).astype(dtype=meta)

    return top_loci, credset_results
