#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
# Modified and updated to process finngen dataset by Miguel Carmona carmona@ebi.ac.uk
# 
#


import sys
import os
import pandas as pd
import json
import subprocess as sp

def main():

    #
    # Args --------------------------------------------------------------------
    #

    study_prefix = 'FINNGEN_R5_'
    # Manifest files from Finngen R5
    in_finngen = 'inputs/r5_finngen.json'
    in_snp_path_list = 'inputs/input_paths_finngen.txt'

    # Path to write main manifest file
    out_manifest = 'finngen.manifest.json'

    # Output directory for individual study fine-mapping results
    output_path = 'output/'

    keep_columns = [
        'code',
        'trait',
        'trait_category',
        'n_cases',
        'n_controls'
    ]

    finngen = (
        pd.read_json(path_or_buf=in_finngen, lines=True)
        .rename(
            columns={
                'phenocode': 'code',
                'phenostring': 'trait',
                'category': 'trait_category',
                'num_cases': 'n_cases',
                'num_controls': 'n_controls'
            }
        )
    )
    finngen = finngen[keep_columns]
    finngen['code'] = study_prefix + finngen['code']
    finngen['n_total'] = finngen['n_cases'] + finngen['n_controls']

    gcs = pd.read_csv(in_snp_path_list, sep='\t', header=None, names=['in_path'])
    gcs['code'] = gcs['in_path'].apply(parse_code, prefix=study_prefix, splitBy='finngen_R5_')

    merged = pd.merge(gcs, finngen, on='code')
    #merged.to_csv('merged.tsv', sep='\t', index=None)

    #
    # Create manifest ---------------------------------------------------------
    #

    with open(out_manifest, 'w') as out_h:
        for in_record in merged.to_dict(orient='records'):
            out_record = {}

            out_record['out_top_loci'] = output_path + in_record['code'] + '/top_loci.json.gz'
            out_record['out_credset'] = output_path + in_record['code'] + '/credible_set.json.gz'
            out_record['out_log'] = output_path + in_record['code'] + '/logfile.txt'

            out_record['in_snp'] = in_record['in_path']
            out_record['study_id'] = in_record['code']
            out_record['n_total'] = int(in_record['n_total'])
            try:
                out_record['n_cases'] = int(in_record['n_cases'])
            except ValueError:
                out_record['n_cases'] = None

            # Write
            out_h.write(json.dumps(out_record) + '\n')

    return 0

def parse_code(in_path, prefix, splitBy):
    ''' Parses the phenotype code from the path name
    '''
    return prefix + in_path.split(splitBy)[-1].split('.')[0]

if __name__ == '__main__':
    main()
