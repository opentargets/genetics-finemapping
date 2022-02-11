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
import gzip

def main():

    # Args
    input_pattern = '/home/ubuntu/results/finemapping/tmp/filtered_input/*.json.gz'
    out_json = 'configs/manifest.json.gz'
    valid_chrom = set([str(chrom) for chrom in range(1, 23)])
    method = 'conditional'

    # Path patterns (local)
    #root = '/Users/em21/Projects/genetics-finemapping'
    # out_path = root + '/output/study_id={0}/phenotype_id={1}/bio_feature={2}/chrom={3}'
    # log_path = root + '/logs/study_id={0}/phenotype_id={1}/bio_feature={2}/chrom={3}'
    # tmp_path = root + '/tmp/study_id={0}/phenotype_id={1}/bio_feature={2}/chrom={3}'
    # ld_ref = '/Users/em21/Projects/reference_data/uk10k_2019Feb/3_liftover_to_GRCh38/output/{chrom}.ALSPAC_TWINSUK.maf01.beagle.csq.shapeit.20131101'
    
    # Path patterns (server)
    # root = '/home/ubuntu/results/finemapping'
    # ld_ref = '/home/ubuntu/data/genotypes/ukb_v3_downsampled10k_plink/ukb_v3_chr{chrom}.downsampled10k'
    root = '/home/js29/genetics-finemapping'
    input_pattern = root + '/tmp/filtered_input/*.json.gz'
    out_path = root + '/output/study_id={0}/phenotype_id={1}/bio_feature={2}/chrom={3}'
    log_path = root + '/logs/study_id={0}/phenotype_id={1}/bio_feature={2}/chrom={3}'
    tmp_path = root + '/tmp/study_id={0}/phenotype_id={1}/bio_feature={2}/chrom={3}'
    
    # In base folder rather than genenetics-finemapping for sharing with coloc pipeline
    #ld_ref = root + '/data/1000Genomes_phase3/EUR/EUR.{chrom}.1000Gp3.20130502'
    ld_ref = '/home/js29/data/ukb_v3_downsampled10k/ukb_v3_chr{chrom}.downsampled10k'

    # Create manifest
    manifest = []
    for in_record in read_json_records(input_pattern):

        # initiate output
        out_record = {}

        # Skip if chromosome is not valid
        if not in_record['chrom'] in valid_chrom:
            continue

        # Add study identifier arguments
        out_record['type'] = in_record.get('type')
        out_record['study_id'] = in_record.get('study_id')
        out_record['phenotype_id'] = in_record.get('phenotype_id', None)
        out_record['bio_feature'] = in_record.get('bio_feature', None)
        out_record['chrom'] = in_record.get('chrom')

        # Add input files
        out_record['in_pq'] = parse_input_name(in_record.get('input_name'))
        out_record['in_ld'] = ld_ref

        # Add output files
        out_record['out_top_loci'] = out_path.format(
            out_record['study_id'], out_record['phenotype_id'],
            out_record['bio_feature'], out_record['chrom']
        ) + '/top_loci.json.gz'
        out_record['out_credset'] = out_path.format(
            out_record['study_id'], out_record['phenotype_id'],
            out_record['bio_feature'], out_record['chrom']
        ) + '/credible_set.json.gz'
        out_record['out_finemap'] = out_path.format(
            out_record['study_id'], out_record['phenotype_id'],
            out_record['bio_feature'], out_record['chrom']
        ) + '/finemap_snp.tsv.gz'
        out_record['out_log'] = log_path.format(
            out_record['study_id'], out_record['phenotype_id'],
            out_record['bio_feature'], out_record['chrom']
        ) + '/logfile.txt'
        out_record['tmpdir'] = tmp_path.format(
            out_record['study_id'], out_record['phenotype_id'],
            out_record['bio_feature'], out_record['chrom']
        )

        # Add method
        out_record['method'] = method
        out_record['pval_threshold'] = in_record.get('pval_threshold')

        manifest.append(out_record)

    # Write manifest as a json
    with gzip.open(out_json, 'w') as out_h:
        for record in manifest:
            out_h.write((json.dumps(record) + '\n').encode())

    return 0

def read_json_records(in_pattern):
    ''' Globs json inputs then yields all records as dicts.
        Expects inputs to be gzipped.
    '''
    for inf in glob(in_pattern):
        with gzip.open(inf, 'r') as in_h:
            for in_record in in_h:
                in_record = json.loads(in_record.decode().rstrip())
                yield in_record

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
