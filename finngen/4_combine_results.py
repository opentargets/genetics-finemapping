#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Jeremy Schwartzentruber
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
    root = '/home/js29/genetics-finemapping/finngen'
    in_top_loci_pattern = root + '/output/*/top_loci.json.gz'
    in_credset_pattern = root + '/output/*/credible_set.json.gz'
    in_logs_pattern = root + '/output/*/logfile.txt'
    out_top_loci = root + '/results/top_loci'
    out_credset = root + '/results/credset'
    out_logs = root + '/results/logfiles.txt'

    # Process top loci
    toploci = spark.read.json(in_top_loci_pattern)
    (
        toploci
        .coalesce(1)
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
    credset = spark.read.json(in_credset_pattern)
    (
        credset
        .repartitionByRange('lead_chrom', 'lead_pos')
        .sortWithinPartitions('lead_chrom', 'lead_pos')
        .write.json(out_credset,
                    compression='gzip',
                    mode='overwrite')
    )
    
    cmd = 'cat {0} | gzip > {1}.gz'.format(in_logs_pattern, out_logs)
    cp = sp.run(cmd, shell=True, stderr=sp.STDOUT)

    return 0

if __name__ == '__main__':

    main()
