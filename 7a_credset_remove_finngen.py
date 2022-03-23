#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Jeremy Schwartzentruber
#
# This code was used for genetics R6 and is no longer needed.
#
import pandas as pd
import pyspark.sql

def main():
    spark = (
        pyspark.sql.SparkSession.builder
        .config("spark.sql.files.ignoreCorruptFiles", "true")
        .config("spark.master", "local[*]")
        .getOrCreate()
    )

    in_credset_pattern = 'finemapping_temp/210923/credset'
    out_credset = 'finemapping_to_merge/210923/credset'
    credset = spark.read.json(in_credset_pattern)
    (
        credset
        .filter(~credset.study_id.contains('FINNGEN'))
        .repartitionByRange('lead_chrom', 'lead_pos')
        .sortWithinPartitions('lead_chrom', 'lead_pos')
        .write.json(out_credset,
                    compression='gzip',
                    mode='overwrite')
    )


if __name__ == '__main__':

    main()
