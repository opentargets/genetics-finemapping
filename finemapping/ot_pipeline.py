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

def run_single_study(in_pq, in_plink,
                     study_id,
                     cell_id=None, gene_id=None, group_id=None, chrom=None,
                     config_file=None,
                     tmp_dir=tempfile.gettempdir(),
                     method='conditional'):
    ''' Runs the top loci and credible set analysis on a single study
    '''

    # Load config file
    with open(config_file, 'r') as in_h:
        config = yaml.load(in_h)
        print('Config: ', config) # DEBUG

    # Load data
    sumstats = finemapping.utils.load_sumstats(
        in_pq,
        study_id,
        cell_id=cell_id,
        gene_id=gene_id,
        group_id=group_id,
        chrom=chrom,
        min_maf=config['min_maf'],
        excl_mhc=config.get('exclude_MHC', None),
        build=config['assembly']
        )

    # Extract top loci
    top_loci = finemapping.top_loci.detect_top_loci(
        sumstats,
        in_plink,
        tmp_dir,
        method=method,
        maf=config['min_maf'],
        cojo_p=float(config['cojo_p']),
        cojo_window=config['cojo_wind'],
        cojo_collinear=config['cojo_colin'],
        clump_dist=config['clump_dist'],
        clump_p=float(config['clump_pval'])
        )

    # Perform credible set analysis on each detected locus
    credset_results = []
    for index_info in top_loci.head(3).to_dict(orient='records'): # DEBUG
        # Run
        credset_res = finemapping.credible_set.run_credible_set_for_locus(
            index_info,
            sumstats,
            top_loci,
            in_plink,
            tmp_dir,
            fm_wind=config['fm_wind'],
            method='conditional'
        )
        # Append result
        credset_results.append(credset_res)

    print(credset_results)

    return 0
