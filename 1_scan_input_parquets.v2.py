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
    spark = (
        pyspark.sql.SparkSession.builder
        .config("spark.master", "local[*]")
        .getOrCreate()
    )
    # sc = spark.sparkContext
    print('Spark version: ', spark.version)

    # Args
    gwas_pval_threshold = 5e-8

    # Paths (local)
    # gwas_pattern = '/Users/em21/Projects/genetics-finemapping/example_data/sumstats/gwas_2/*.parquet'
    # mol_pattern = '/Users/em21/Projects/genetics-finemapping/example_data/sumstats/molecular_trait_2/*.parquet'

    # Paths (server)
    gwas_pattern = '/home/ubuntu/data/sumstats/filtered/significant_window_2mb/gwas/*.parquet'
    mol_pattern = '/home/ubuntu/data/sumstats/filtered/significant_window_2mb/molecular_trait/*.parquet'
    out_path = '/home/ubuntu/results/finemapping/tmp/filtered_input'

    # Stop if path exists
    if os.path.exists(out_path):
        sys.exit('Error: output path exists {}'.format(out_path))

    #
    # Process GWAS
    #

    for inf in glob(gwas_pattern):
        
        # Load
        inf = os.path.abspath(inf)
        df = (
            spark.read.parquet(inf)
                .withColumn('pval_threshold', lit(gwas_pval_threshold))
                .withColumn('input_name', lit(inf))
        )

        # Process
        df = (
            df.filter(col('pval') < col('pval_threshold'))
            .select('type', 'study_id', 'phenotype_id', 'bio_feature', 'gene_id', 'chrom', 'pval_threshold', 'input_name')
            .distinct()
        )

        # Save
        (
            df
            .coalesce(1)
            .write.json(out_path,
                        compression='gzip',
                        mode='append')
        )

        df.unpersist()
    
    #
    # Process molecular traits
    #

    for inf in glob(mol_pattern):

        # Load
        inf = os.path.abspath(inf)
        df = (
            spark.read.parquet(inf)
            .withColumn('pval_threshold', (0.05 / col('num_tests')))
            .withColumn('pval_threshold', when(col('pval_threshold') > gwas_pval_threshold,
                                            col('pval_threshold'))
                        .otherwise(gwas_pval_threshold))
            .drop('num_tests')
            .withColumn('input_name', lit(inf))
        )
        
        # Process
        df = (
            df.filter(col('pval') < col('pval_threshold'))
            .select('type', 'study_id', 'phenotype_id', 'bio_feature', 'gene_id', 'chrom', 'pval_threshold', 'input_name')
            .distinct()
        )

        # Save
        (
            df
            .coalesce(1)
            .write.json(out_path,
                        compression='gzip',
                        mode='append')
        )

    return 0

if __name__ == '__main__':

    main()
