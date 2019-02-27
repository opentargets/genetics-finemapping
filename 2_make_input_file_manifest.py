#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import json
import os
from pprint import pprint
import pandas as pd
from numpy import nan
from glob import glob

def main():

    # Args
    input_json = glob('tmp/filtered_input.json/*.json')[0]
    out_manifest = 'configs/input_files.config.tsv'
    # valid_chrom = set([str(chrom) for chrom in range(1, 23)])
    valid_chrom = set(['22'])
    method = 'conditional'

    # Path patterns
    out_path = '/Users/em21/Projects/genetics-finemapping/output/study_id={0}/phenotype_id={1}/biofeature={2}/chrom={3}'
    log_path = '/Users/em21/Projects/genetics-finemapping/logs/study_id={0}/phenotype_id={1}/biofeature={2}/chrom={3}'
    tmp_path = '/Users/em21/Projects/genetics-finemapping/tmp/study_id={0}/phenotype_id={1}/biofeature={2}/chrom={3}'
    ld_ref = '/Users/em21/Projects/reference_data/uk10k_2019Feb/3_liftover_to_GRCh38/output/{chrom}.ALSPAC_TWINSUK.maf01.beagle.csq.shapeit.20131101'
    
    # Create manifest
    manifest = []
    with open(input_json, 'r') as in_h:
        for in_record in in_h:
            in_record = json.loads(in_record)
            out_record = {}

            # Add study identifier arguments
            out_record['type'] = in_record.get('type')
            out_record['study_id'] = in_record.get('study_id')
            out_record['phenotype_id'] = in_record.get('phenotype_id', None)
            out_record['biofeature'] = in_record.get('biofeature', None)
            out_record['chrom'] = in_record.get('chrom')

            # Add input files
            out_record['in_pq'] = parse_input_name(in_record.get('input_name'))
            out_record['in_ld'] = ld_ref

            # Add output files
            out_record['out_top_loci'] = out_path.format(
                out_record['study_id'], out_record['phenotype_id'],
                out_record['biofeature'], out_record['chrom']
            ) + '/top_loci.json.gz'
            out_record['out_credset'] = out_path.format(
                out_record['study_id'], out_record['phenotype_id'],
                out_record['biofeature'], out_record['chrom']
            ) + '/credible_set.json.gz'
            out_record['out_log'] = log_path.format(
                out_record['study_id'], out_record['phenotype_id'],
                out_record['biofeature'], out_record['chrom']
            ) + '/logfile.txt'
            out_record['tmpdir'] = tmp_path.format(
                out_record['study_id'], out_record['phenotype_id'],
                out_record['biofeature'], out_record['chrom']
            )

            # Add method
            out_record['method'] = method
            out_record['pval_threshold'] = in_record.get('pval_threshold')

            manifest.append(out_record)

    # Convert to a dataframe
    df = pd.DataFrame(manifest)
    col_order = ['type', 'in_pq', 'in_ld', 'study_id', 'phenotype_id',
                 'biofeature', 'chrom', 'out_top_loci', 'out_credset',
                 'out_log', 'tmpdir', 'method', 'pval_threshold']
    df = df.loc[:, col_order]

    # Write
    df.to_csv(out_manifest, sep='\t', header=None, index=None, na_rep='None')

    return 0

def parse_input_name(s):
    ''' Parses the required input name. Spark's input_file_name() returns the
        nested parquet file, I need the top level parquet.
    '''
    # Strip nested parquet
    out_s = s.split('.parquet')[0] + '.parquet'
    # Stip file://
    out_s = out_s.replace('file://', '')
    return out_s

if __name__ == '__main__':

    main()
