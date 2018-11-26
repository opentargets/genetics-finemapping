#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import json
import os
import dask.dataframe as dd
from pprint import pprint
import pandas as pd
from numpy import nan

def main():

    # Args
    data_pattern = 'input/{type}/*/*.parquet'
    cols = ['study_id', 'trait_id', 'group_id', 'cell_id', 'chrom']
    valid_chrom = set([str(chrom) for chrom in range(1, 23)])
    out_manifest = 'configs/input_files.config.tsv'

    manifest = []

    # for type in ['gwas']::
    for type in ['gwas', 'molecular_qtl']:
        # Load list of all studies
        in_path = data_pattern.format(type=type)
        studies = dd.read_parquet(
            path=in_path,
            columns=cols,
            engine='fastparquet'
        )
        # Drop duplicates
        studies = studies.drop_duplicates()
        # Convert to pandas
        studies = studies.compute()
        # Add manifest record
        for i, row in studies.iterrows():
            # Skip if not valid chrom
            if not row.chrom in valid_chrom:
                continue
            # Populate with study specific info
            record = {}
            for key in cols:
                record[key] = row[key] if row[key] else ''
            # Populate other fields
            record['in_pq'] = os.path.abspath('input/{0}/{1}'.format(type, record['study_id']))
            record['out_top_loci'] = os.path.abspath('output/study_id={0}/cell_id={1}/group_id={2}/trait_id={3}/chrom={4}/top_loci.parquet'.format(record['study_id'], record['cell_id'], record['group_id'], record['trait_id'], record['chrom']))
            record['out_credset'] = os.path.abspath('output/study_id={0}/cell_id={1}/group_id={2}/trait_id={3}/chrom={4}/credible_set.parquet'.format(record['study_id'], record['cell_id'], record['group_id'], record['trait_id'], record['chrom']))
            record['out_log'] = os.path.abspath('logs/study_id={0}/cell_id={1}/group_id={2}/trait_id={3}/chrom={4}/logfile.txt'.format(record['study_id'], record['cell_id'], record['group_id'], record['trait_id'], record['chrom']))
            record['tmpdir'] = os.path.abspath('tmp/study_id={0}/cell_id={1}/group_id={2}/trait_id={3}/chrom={4}/'.format(record['study_id'], record['cell_id'], record['group_id'], record['trait_id'], record['chrom']))
            record['ld_ref'] = os.path.abspath('input/ld/EUR.{chrom}.1000Gp3.20130502')
            record['method'] = 'conditional'
            # Add manifest
            manifest.append(record)

    # Convert record to a dataframe
    df = pd.DataFrame(manifest)
    col_order = ['in_pq', 'ld_ref', 'study_id', 'cell_id', 'group_id', 'trait_id', 'chrom', 'method', 'out_top_loci', 'out_credset', 'out_log', 'tmpdir']
    df = df.loc[:, col_order]
    df = df.replace('', nan)

    # Write
    df.to_csv(out_manifest, sep='\t', header=None, index=None, na_rep='None')

    return 0

if __name__ == '__main__':

    main()
