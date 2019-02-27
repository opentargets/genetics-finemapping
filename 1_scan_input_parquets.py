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

def main():

    # Args
    gwas_pval_threshold = 5e-8

    # Make spark session
    spark = pyspark.sql.SparkSession.builder.getOrCreate()
    # sc = spark.sparkContext
    print('Spark version: ', spark.version)

    # Load GWAS and molecular trait separately, then union
    gwas = (
        spark.read.parquet('example_data/sumstats/gwas/*.parquet')
             .withColumn('pval_threshold', lit(gwas_pval_threshold))
             .withColumn('input_name', input_file_name())
    )
    mol = (
        spark.read.parquet('example_data/sumstats/molecular_trait/*.parquet')
            .withColumn('pval_threshold', (0.05 / col('num_tests')))
            .withColumn('pval_threshold', when(col('pval_threshold') > gwas_pval_threshold,
                                               col('pval_threshold'))
                                          .otherwise(gwas_pval_threshold))
             .drop('num_tests')
             .withColumn('input_name', input_file_name())
    )
    df = gwas.unionByName(mol)
    
    # Process
    df = (
        df.filter(col('pval') < col('pval_threshold'))
          .select('type', 'study_id', 'phenotype_id', 'biofeature', 'gene_id', 'chrom', 'pval_threshold', 'input_name')
          .filter(col('chrom') == '22') # DEBUG
          .distinct()
    )

    # Write
    (
        df.coalesce(1)
          .write.json('tmp/filtered_input.json',
                      mode='overwrite')
    )

    # input_file_name

    return 0

if __name__ == '__main__':

    main()
