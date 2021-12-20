#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Jeremy Schwartzentruber
#
# This script was used when we found that some top_loci were duplicated, which was
# due to eQTL catalogue ingestion before duplicate SNPs were removed. In general
# this script should not be needed.
#
'''
# Set SPARK_HOME and PYTHONPATH to use 2.4.0
export PYSPARK_SUBMIT_ARGS="--driver-memory 8g pyspark-shell"
export SPARK_HOME=/Users/em21/software/spark-2.4.0-bin-hadoop2.7
export PYTHONPATH=$SPARK_HOME/python:$SPARK_HOME/python/lib/py4j-2.4.0-src.zip:$PYTHONPATH
'''

import pyspark.sql
from pyspark.sql.types import *
from pyspark.sql.functions import *
import os
from shutil import copyfile
from glob import glob
import yaml
import subprocess as sp

def main():

    # Make spark session
    # Using `ignoreCorruptFiles` will skip empty files
    spark = (
        pyspark.sql.SparkSession.builder
        .config("spark.sql.files.ignoreCorruptFiles", "true")
        .config("spark.master", "local[*]")
        .getOrCreate()
    )
    # sc = spark.sparkContext
    print('Spark version: ', spark.version)

    # Args
    #root = '/home/ubuntu/results/finemapping'
    root = '/home/js29/genetics-finemapping'
    #in_top_loci_pattern = root + '/output/study_id=*/phenotype_id=*/bio_feature=*/chrom=*/top_loci.json.gz'
    #in_credset_pattern = root + '/output/study_id=*/phenotype_id=*/bio_feature=*/chrom=*/credible_set.json.gz'
    in_top_loci_pattern = root + '/top_loci.dedup.json.gz'
    in_credset_pattern = root + '/credible_set.dedup.json.gz'
    
    out_top_loci = root + '/results/top_loci'
    out_credset = root + '/results/credset'

    with open('configs/analysis.config.yaml', 'r') as in_h:
        config_dict = yaml.load(in_h, Loader=yaml.FullLoader)

    # Process top loci
    toploci = spark.read.json(in_top_loci_pattern)
    # Cols phenotype_id and bio_feature might not be present if no molecular
    # QTL data were processed
    if not 'phenotype_id' in toploci.columns:
        toploci = toploci.withColumn('phenotype_id', lit(None))
    if not 'bio_feature' in toploci.columns:
        toploci = toploci.withColumn('bio_feature', lit(None))

    # Remove any duplicate rows of top_loci
    last_n = toploci.count()
    toploci = toploci.dropDuplicates(subset=['study_id', 'bio_feature', 'phenotype_id', 'type', 'variant_id'])
    num_toploci = toploci.count()
    print('{} toploci removed that were duplicates'.format( last_n - num_toploci ))

    num_gwas_loci = toploci.filter(col('type') == 'gwas').count()
    num_qtl_loci = num_toploci - num_gwas_loci
    print(f'{num_toploci} toploci remain. {num_gwas_loci} GWAS and {num_qtl_loci} QTL')

    (
        toploci
        .coalesce(1)
        .orderBy('study_id', 'phenotype_id', 'bio_feature',
                 'chrom', 'pos')
        .write.json(out_top_loci,
                    compression='gzip',
                    mode='overwrite')
    )
    
    # Copy to single file (this seems to only copy the first file in the folder...)
    copyfile(
        glob(out_top_loci + '/part-*.json.gz')[0],
        out_top_loci + '.json.gz'
    )
    
    # Process cred set
    credset = spark.read.json(in_credset_pattern)
    if not 'phenotype_id' in credset.columns:
        credset = credset.withColumn('phenotype_id', lit(None))
    if not 'bio_feature' in credset.columns:
        credset = credset.withColumn('bio_feature', lit(None))

    # Remove any duplicate rows
    last_n = credset.count()
    credset = credset.dropDuplicates(subset=['study_id', 'bio_feature', 'phenotype_id', 'lead_variant_id', 'tag_variant_id'])
    num_credset_rows = credset.count()
    print('{} credset rows removed that were duplicates'.format( last_n - num_credset_rows ))
    print(f'{num_credset_rows} credset rows remain.')

    (
        credset
        .repartitionByRange('lead_chrom', 'lead_pos')
        .sortWithinPartitions('lead_chrom', 'lead_pos')
        .write.json(out_credset,
                    compression='gzip',
                    mode='overwrite')
    )
    
    if config_dict['run_finemap']:
        # Process finemap output
        in_finemap_pattern = root + '/output/study_id=*/phenotype_id=*/bio_feature=*/chrom=*/finemap_snp.tsv.gz'
        out_finemap_snp = root + '/results/finemap_snp'
        finemap_res = (
            spark.read.option("delimiter", "\t")
            .option("header", "true")
            .csv(in_finemap_pattern)
        )
        (
            finemap_res
            .filter(col('prob') > 0.0001)
            .repartitionByRange('chromosome', 'position')
            .sortWithinPartitions('chromosome', 'position')
            .write.csv(out_finemap_snp + "_filtered",
                    header=True,
                    compression='gzip',
                    mode='overwrite')
        )
        (
            finemap_res
            .repartitionByRange('chromosome', 'position')
            .sortWithinPartitions('chromosome', 'position')
            .write.csv(out_finemap_snp,
                    header=True,
                    compression='gzip',
                    mode='overwrite')
        )

        cmd = 'zcat {0}/part-00000*.csv.gz | head -n 1 > {0}.csv'.format(out_finemap_snp)
        cp = sp.run(cmd, shell=True, stderr=sp.STDOUT)

        cmd = 'for f in {0}/part-*.csv.gz; do zcat $f | sed "1d" >> {0}.csv; done; gzip -f {0}.csv'.format(out_finemap_snp)
        cp = sp.run(cmd, shell=True, stderr=sp.STDOUT)

    return 0

if __name__ == '__main__':

    main()
