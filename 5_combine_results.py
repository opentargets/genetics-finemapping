#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
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
import yaml
import os
from shutil import copyfile
from glob import glob
import yaml
import subprocess as sp

def main():
    # Load analysis config file
    config_file = 'configs/analysis.config.yaml'
    with open(config_file, 'r') as in_h:
        config_dict = yaml.safe_load(in_h)

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
    in_top_loci_pattern = os.path.join(config_dict['finemapping_output_dir'],
                                       'output/study_id=*/phenotype_id=*/bio_feature=*/chrom=*/top_loci.json.gz')
    in_credset_pattern = os.path.join(config_dict['finemapping_output_dir'],
                                      'output/study_id=*/phenotype_id=*/bio_feature=*/chrom=*/credible_set.json.gz')

    out_top_loci = os.path.join(config_dict['finemapping_output_dir'], 'results/top_loci')
    if not os.path.exists(out_top_loci):
        os.makedirs(out_top_loci)

    out_credset = os.path.join(config_dict['finemapping_output_dir'], 'results/credset')
    if not os.path.exists(out_credset):
        os.makedirs(out_credset)

    with open('configs/analysis.config.yaml', 'r') as in_h:
        config_dict = yaml.safe_load(in_h)

    # Process top loci
    toploci = spark.read.json(in_top_loci_pattern)
    # Cols phenotype_id and bio_feature might not be present if no molecular
    # QTL data were processed
    if not 'phenotype_id' in toploci.columns:
        toploci = toploci.withColumn('phenotype_id', lit(None))
    if not 'bio_feature' in toploci.columns:
        toploci = toploci.withColumn('bio_feature', lit(None))
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
        in_finemap_pattern = os.path.join(config_dict['finemapping_output_dir'], '/output/study_id=*/phenotype_id=*/bio_feature=*/chrom=*/finemap_snp.tsv.gz')
        out_finemap_snp = os.path.join(config_dict['finemapping_output_dir'], '/results/finemap_snp')
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
