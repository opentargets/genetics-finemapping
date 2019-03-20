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
from functools import reduce

def main():

    # Make spark session
    spark = pyspark.sql.SparkSession.builder.getOrCreate()
    # sc = spark.sparkContext
    print('Spark version: ', spark.version)

    # Args
    gwas_pval_threshold = 0.05

    # Paths (local)
    # gwas_pattern = '/Users/em21/Projects/genetics-finemapping/example_data/sumstats/gwas_2/*.parquet'
    # mol_pattern = '/Users/em21/Projects/genetics-finemapping/example_data/sumstats/molecular_trait_2/*.parquet'
    # out_path = '/Users/em21/Projects/genetics-finemapping/example_data/sumstats_filtered'

    # Paths (server)
    gwas_pattern = '/home/em21/data/sumstats/gwas_2/*.parquet'
    mol_pattern = '/home/em21/data/sumstats/molecular_trait_2/*.parquet'
    out_path = '/home/em21/data/sumstats_filtered'

    # Load GWAS dfs
    gwas_dfs = []
    for inf in glob(gwas_pattern):
        inf = os.path.abspath(inf)
        df = spark.read.parquet(inf)
        gwas_dfs.append(df)
    
    # Load molecular trait dfs
    mol_dfs = []
    for inf in glob(mol_pattern):
        inf = os.path.abspath(inf)
        df = (
            spark.read.parquet(inf)
            .drop('num_tests')
        )
        mol_dfs.append(df)

    # Take union
    df = reduce(
        pyspark.sql.DataFrame.unionByName,
        gwas_dfs + mol_dfs
    )
    
    # Process
    df = (
        df.filter(col('pval') < gwas_pval_threshold)
          .repartitionByRange('chrom', 'pos')
    )

    # Write
    (
        df.write.json(out_path,
                      compression='gzip',
                      mode='overwrite')
    )

    # input_file_name

    return 0

if __name__ == '__main__':

    main()
