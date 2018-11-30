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
                     cell_id=None,
                     group_id=None,
                     trait_id=None,
                     chrom=None,
                     analysis_config=None,
                     tmp_dir=tempfile.gettempdir(),
                     method='conditional',
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
        cell_id=cell_id,
        group_id=group_id,
        trait_id=trait_id,
        chrom=chrom,
        min_maf=analysis_config['min_maf'],
        excl_mhc=analysis_config.get('exclude_MHC', None),
        build=analysis_config['assembly'],
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
        cojo_p=float(analysis_config['cojo_p']),
        cojo_window=analysis_config['cojo_wind'],
        cojo_collinear=analysis_config['cojo_colin'],
        clump_dist=analysis_config['clump_dist'],
        clump_p=float(analysis_config['clump_pval']),
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
    for c,d in zip(columns, dtypes):
        df[c] = pd.Series(dtype=d)
    return df

# def run_all_studies(gwas_sumstats_dir, mol_sumstats_dir, manifest_file,
#         results_dir, analysis_config_file, temp_dir, skip_existing,
#         existing_results_dir, n_cores=1):
#     ''' Run the fine-mapping pipeline on all studies in the manifest file.
#         Optionally, skips studies that already have results.
#     Args:
#         gwas_sumstats_dir (dir): directory containing gwas summary stats
#         mol_sumstats_dir (dir): dir containing molecular trait gwas sumstats
#         manifest_file (file): manifest file specifying which studies to run
#         plink_ld (str): pattern matching ld reference for european samples
#         results_dir (dir): where to save the results
#         analysis_config_file (file): config file containing analysis parameters
#         temp_dir (dir): where to store temporary files
#         skip_completed (bool): whether studies with existing results should be
#                                skipped
#         skipped_results (dir): location of results for skipped studies
#         n_cores (int): number of cores to use
#     '''
#
#     #
#     # Setup
#     #
#
#     # Set scheduler params
#     client = Client(n_workers=n_cores)
#
#     # Load the manifest
#     with open(manifest_file, 'r') as in_h:
#         manifest = json.load(in_h)
#
#     # DEBUG
#     print('Input manifest: ', end='')
#     pprint(manifest)
#
#     # TEMP solution, output directory must be local
#     assert finemapping.utils.is_local(results_dir)
#
#     #
#     # Conpare input manifest to existing results if skip_existing==True
#     #
#
#     if skip_existing is True:
#
#         print('Skipping existing...')
#
#         # Load existing results files
#         existing_top_loci = dask.dataframe.read_parquet(
#             os.path.join(existing_results_dir, 'top_loci.parquet'),
#             engine='fastparquet').compute()
#         existing_cred_set = dask.dataframe.read_parquet(
#             os.path.join(existing_results_dir, 'credible_sets.parquet'),
#             engine='fastparquet').compute()
#
#         # Load existing
#         existing_manifest_file = os.path.join(
#             existing_results_dir, 'manifest.completed.json')
#         with open(existing_manifest_file, 'r') as in_h:
#             existing_manifest = json.load(in_h)
#
#         # Compare manifest to existing files
#         manifest = compare_existing_with_manifest(manifest, existing_manifest)
#
#
#     # TEMP solution. Check that there is at least 1 study to run.
#     assert len(manifest) >= 1
#
#     #
#     # Run the pipeline
#     #
#
#     # Create delayed function
#     run_single_study_delayed = dask.delayed(run_single_study, nout=2)
#
#     # Initiate empty results lists
#     top_loci_results_list = []
#     cred_set_results_list = []
#
#     # Created delayed execution for each study in the manfiest
#     for study_info in manifest:
#
#         # Create pq input variable
#         if study_info['type'] == 'gwas':
#             in_pq = os.path.join(gwas_sumstats_dir, study_info['study_id'])
#         elif study_info['type'] == 'molecular_qtl':
#             in_pq = os.path.join(mol_sumstats_dir, study_info['study_id'])
#
#         # Create temp input variable
#         tmp_dir = os.path.join(temp_dir, study_info['study_id'])
#
#         # Submit delayed
#         top_loci, credset = run_single_study_delayed(
#             in_pq=in_pq,
#             in_plink=study_info['ld_ref'],
#             study_id=study_info['study_id'],
#             cell_id=study_info['cell_id'],
#             group_id=study_info['group_id'],
#             trait_id=study_info['trait_id'],
#             chrom=study_info['chrom'],
#             analysis_config_file=analysis_config_file,
#             tmp_dir=tmp_dir,
#             method=study_info['method'] )
#
#         # Add results to lists
#         top_loci_results_list.append(top_loci)
#         cred_set_results_list.append(credset)
#
#     # Make dask.dfs from the delayed objects
#     top_loci_dd = dask.dataframe.from_delayed(
#         top_loci_results_list,
#         meta=finemapping.utils.get_meta_info(type='top_loci'))
#     cred_set_dd = dask.dataframe.from_delayed(
#         cred_set_results_list,
#         meta=finemapping.utils.get_meta_info(type='cred_set'))
#
#     # Compute the results. Brings everything into memory on single machine.
#     print('Starting compute')
#     top_loci_df, cred_set_df = dask.compute(top_loci_dd, cred_set_dd, scheduler='single-threaded') # DEBUG
#
#     # Concatenate existing results to new results
#     if skip_existing is True:
#         top_loci_df = pd.concat([top_loci_df, existing_top_loci],
#                                 ignore_index=True)
#         cred_set_df = pd.concat([cred_set_df, existing_cred_set],
#                                 ignore_index=True)
#
#     #
#     # Write results
#     #
#
#     # Make output folder
#     os.makedirs(results_dir, exist_ok=True)
#
#     # Write
#     top_loci_df.to_parquet(
#         fname=os.path.join(results_dir, 'top_loci.parquet'),
#         engine='fastparquet',
#         compression='snappy'
#     )
#     cred_set_df.to_parquet(
#         fname=os.path.join(results_dir, 'credible_sets.parquet'),
#         engine='fastparquet',
#         compression='snappy'
#     )
#
#     #
#     # Write completed manifest
#     #
#
#     # Add todays date to manifest records
#     dt = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
#     for record in manifest:
#         record['completed_date'] = dt
#
#     # Merge existing with current manifest
#     if skip_existing is True:
#         manifest = manifest + existing_manifest
#
#     # Write manifest
#     out_manifest_file = os.path.join(
#         results_dir, 'manifest.completed.json')
#     with open(out_manifest_file, 'w') as out_h:
#         json.dump(manifest, out_h, indent=2)
#
#     return 0

# def compare_existing_with_manifest(manifest, existing_manifest):
#     ''' Removes existing records from the manifest
#     Args:
#         manifest (list of dicts)
#         existing_manifest (list of dicts)
#     Returns:
#         manifest with existing studies removed (list of dicts)
#     '''
#
#     # Make set of existing keys
#     existing_keys = set(
#         [make_manifest_record_key(record) for record in existing_manifest]
#     )
#
#     # Make list of non-existing records
#     novel_manifest = []
#     for record in manifest:
#         if make_manifest_record_key(record) not in existing_keys:
#             novel_manifest.append(record)
#
#     return novel_manifest
#
# def make_manifest_record_key(record):
#     ''' Given a single record from the manifest, return unique key
#     '''
#     fields = ['study_id', 'trait_id', 'cell_id', 'group_id', 'chrom']
#     return '-'.join([record[field] if record[field] else '' for field in fields])
