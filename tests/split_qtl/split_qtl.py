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
from glob import glob

def main():

    # Args
    mol_pattern = '/home/emountjoy_statgen/data/sumstats/molecular_trait/*.parquet'
    out_dir = '/home/emountjoy_statgen/data/sumstats/molecular_trait_split/'

    # Make spark session
    spark = pyspark.sql.SparkSession.builder.getOrCreate()
    print('Spark version: ', spark.version)

    # Process each
    for inf in glob(mol_pattern):

        # Load
        df = spark.read.parquet(inf)

        # Write
        outf = os.path.join(out_dir, os.path.basename(inf))
        (
            df.write
            .partitionBy('biofeature', 'chrom')
            .parquet(
                outf,
                mode='overwrite',
                compression='snappy'
            )
        )
    
    return 0

if __name__ == '__main__':

    main()
