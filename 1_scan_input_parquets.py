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

    # Paths
    gwas_pattern = '/home/js29/genetics-finemapping/data/filtered/significant_window_2mb/gwas/*.parquet'
    mol_pattern = '/home/js29/genetics-finemapping/data/filtered/significant_window_2mb/molecular_trait/*.parquet'
    out_path = '/home/js29/genetics-finemapping/tmp/filtered_input'

    # Load GWAS dfs
    strip_path_gwas = udf(lambda x: x.replace('file:', '').split('/part-')[0], StringType())
    gwas_dfs = (
        spark.read.parquet(gwas_pattern)
            .withColumn('pval_threshold', lit(gwas_pval_threshold))
            .withColumn('input_name', strip_path_gwas(input_file_name()))
    )
    
    # Load molecular trait dfs
    # This has to be done separately, followed by unionByName as the hive
    # parititions differ across datasets due to different tissues
    # (bio_features) and chromosomes
    strip_path_mol = udf(lambda x: x.replace('file:', ''), StringType())
    mol_dfs = []
    for inf in glob(mol_pattern):
        df = (
            spark.read.parquet(inf)
            .withColumn('pval_threshold', (0.05 / col('num_tests')))
            .withColumn('pval_threshold', when(col('pval_threshold') > gwas_pval_threshold,
                                            col('pval_threshold'))
                        .otherwise(gwas_pval_threshold))
            .drop('num_tests')
            .withColumn('input_name', strip_path_mol(lit(inf)))
        )
        mol_dfs.append(df)

    # Take union
    df = reduce(
        pyspark.sql.DataFrame.unionByName,
        [gwas_dfs] + mol_dfs
    )
    
    # Process
    df = (
        df.filter(col('pval') < col('pval_threshold'))
          .select('type', 'study_id', 'phenotype_id', 'bio_feature', 'gene_id', 'chrom', 'pval_threshold', 'input_name')
          .distinct()
    )

    # Write
    (
        df
          .coalesce(300)
          .write.json(out_path,
                      compression='gzip',
                      mode='overwrite')
    )

    return 0

if __name__ == '__main__':

    main()
