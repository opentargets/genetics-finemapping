#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Jeremy Schwartzentruber
#
import pandas as pd
import pyspark.sql

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

    spark = (
        pyspark.sql.SparkSession.builder
        .config("spark.sql.files.ignoreCorruptFiles", "true")
        .config("spark.master", "local[*]")
        .getOrCreate()
    )

    in_credset_pattern = 'finemapping_results/*/credset'
    out_credset = 'finemapping_merged/credset'
    credset = spark.read.json(in_credset_pattern)
    (
        credset
        .repartitionByRange('lead_chrom', 'lead_pos')
        .sortWithinPartitions('lead_chrom', 'lead_pos')
        .write.json(out_credset,
                    compression='gzip',
                    mode='overwrite')
    )


if __name__ == '__main__':

    main()
