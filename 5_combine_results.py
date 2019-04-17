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
import os
from shutil import copyfile
from glob import glob

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
    in_top_loci_pattern = 'output/study_id=*/phenotype_id=*/bio_feature=*/chrom=*/top_loci.json.gz'
    in_credset_pattern = 'output/study_id=*/phenotype_id=*/bio_feature=*/chrom=*/credible_set.json.gz'
    out_top_loci = 'results/top_loci'
    out_credset = 'results/credset'

    # Process top loci 
    (
        spark.read.json(in_top_loci_pattern)
        .coalesce(1)
        .orderBy('study_id', 'phenotype_id', 'bio_feature',
                 'chrom', 'pos')
        .write.json(out_top_loci,
                    compression='gzip',
                    mode='overwrite')
    )
    
    # Copy to single file
    copyfile(
        glob(out_top_loci + '/part-*.json.gz')[0],
        out_top_loci + '.json.gz'
    )
    
    # Process cred set
    (
        spark.read.json(in_credset_pattern)
        .repartitionByRange('study_id', 'phenotype_id', 'bio_feature',
                            'lead_chrom', 'lead_pos')
        .write.json(out_credset,
                    compression='gzip',
                    mode='overwrite')
    )
    


    return 0

if __name__ == '__main__':

    main()
