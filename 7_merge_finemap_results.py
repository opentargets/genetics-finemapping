#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Jeremy Schwartzentruber
#
import os
import argparse
import pandas as pd
import pyspark.sql
from pyspark.sql.functions import *
import subprocess as sp

def main():
    # Don't merge top_loci in this way, because it causes all loci to have the same
    # fields - i.e. GWAS top loci have "bio_feature" fields, which become null, and
    # then the JSON can't be read in.
    #top_loci0 = pd.read_json('finemapping_results/190612/top_loci.json.gz', orient='records', lines=True)
    #top_loci1 = pd.read_json('finemapping_results/210309/top_loci.json.gz', orient='records', lines=True)
    #top_loci2 = pd.read_json('finemapping_results/210421/top_loci.json.gz', orient='records', lines=True)
    #top_loci3 = pd.read_json('finemapping_results/finngen_210509/top_loci.json.gz', orient='records', lines=True)
    #top_loci4 = pd.read_json('finemapping_results/210507/top_loci.json.gz', orient='records', lines=True)

    #all_top_loci = pd.concat([top_loci0, top_loci1, top_loci2, top_loci3])

    #all_top_loci.to_json(
    #    'finemapping_results/top_loci.json.gz',
    #    orient='records',
    #    lines=True,
    #    compression='gzip',
    #    double_precision=15
    #)
    args = parse_args()

    spark = (
        pyspark.sql.SparkSession.builder
        .config("spark.sql.files.ignoreCorruptFiles", "true")
        .config("spark.master", "local[*]")
        .getOrCreate()
    )

    # Read in the top loci
    top_loci_prev = spark.read.json(os.path.join(args.prev_results, 'top_loci.json.gz'))
    top_loci_new = spark.read.json(os.path.join(args.new_results, 'top_loci.json.gz'))

    nrows_prev_start = top_loci_prev.count()
    print(f'Rows in previous top_loci: {nrows_prev_start}')
    print('Rows in new top_loci: {0}'.format(top_loci_new.count()))

    if args.remove_previous_finngen:
        # Remove FinnGen rows from the previous results
        top_loci_prev = top_loci_prev.filter(~col('study_id').contains('FINNGEN'))
        nrows_prev_end = top_loci_prev.count()
        print('FinnGen rows removed from previous top_loci: {0}'.format(nrows_prev_start - nrows_prev_end))

    top_loci = top_loci_prev.unionByName(top_loci_new, allowMissingColumns=True)

    # Write out the top loci
    print('Writing top_loci rows: {0}'.format(top_loci.count()))
    top_loci_folder = os.path.join(args.output, 'top_loci')
    (
        top_loci
        .write.json(top_loci_folder, 
                    compression='gzip',
                    mode='overwrite')
    )
    cmd = 'zcat {0}/part*.json.gz | gzip > {0}.json.gz'.format(top_loci_folder)
    cp = sp.run(cmd, shell=True, stderr=sp.STDOUT)

    # copyfile(
    #     glob(os.path.join(args.output, 'top_loci') + '/part-*.json.gz')[0],
    #     os.path.join(args.output, 'top_loci.json.gz')
    # )

    # Read in the credsets
    credset_prev = spark.read.json(os.path.join(args.prev_results, 'credset'))
    credset_new = spark.read.json(os.path.join(args.new_results, 'credset'))

    nrows_prev_start = credset_prev.count()
    print(f'Rows in previous credsets: {nrows_prev_start}')
    print('Rows in new credsets: {0}'.format(credset_new.count()))

    if args.remove_previous_finngen:
        # Remove FinnGen rows from the previous results
        credset_prev = credset_prev.filter(~col('study_id').contains('FINNGEN'))
        nrows_prev_end = credset_prev.count()
        print('FinnGen rows removed from previous credsets: {0}'.format(nrows_prev_start - nrows_prev_end))

    credset = credset_prev.unionByName(credset_new, allowMissingColumns=True)

    # Write out the credsets
    print('Writing credset rows: {0}'.format(credset.count()))
    (
        credset
        .repartitionByRange('lead_chrom', 'lead_pos')
        .sortWithinPartitions('lead_chrom', 'lead_pos')
        .write.json(os.path.join(args.output, 'credset'),
                    compression='gzip',
                    mode='overwrite')
    )

    # in_credset_pattern = 'finemapping_to_merge/*/credset'
    # out_credset = 'finemapping_merged/credset'
    # credset = spark.read.json(in_credset_pattern)
    # (
    #     credset
    #     .repartitionByRange('lead_chrom', 'lead_pos')
    #     .sortWithinPartitions('lead_chrom', 'lead_pos')
    #     .write.json(out_credset,
    #                 compression='gzip',
    #                 mode='overwrite')
    # )


def parse_args():
    ''' Load command line args
    '''
    p = argparse.ArgumentParser()

    # Add input files
    p.add_argument('--prev_results',
                   metavar="<file>",
                   help=("Input: previous fine-mapping results folder"),
                   type=str,
                   required=True)

    p.add_argument('--new_results',
                   metavar="<file>",
                   help=("Input: new fine-mapping results folder"),
                   type=str,
                   required=True)

    p.add_argument('--output',
                   metavar="<file>",
                   help=("Output root folder"),
                   type=str,
                   required=True)
    
    p.add_argument('--remove_previous_finngen',
                   help=("If set, remove FinnGen rows from the previous results"),
                   action='store_true')

    args = p.parse_args()

    return args


if __name__ == '__main__':

    main()
