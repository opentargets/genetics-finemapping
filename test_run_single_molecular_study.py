#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import finemapping.ot_pipeline
import os

def main():


    # Set up test params
    tmp_dir = 'tmp/GTEX7'
    out_dir = 'output/GTEX7'

    # Create output tmp dir
    # os.makedirs(tmp_dir, exist_ok=True)

    # Run example
    top_loci, credsets = finemapping.ot_pipeline.run_single_study(
        in_pq='input/molecular_qtl/GTEX7',
        in_plink='input/ld/EUR.{chrom}.1000Gp3.20130502',
        study_id='GTEX7',
        analysis_config_file='configs/analysis.config.yaml',
        chrom='1',
        cell_id='UBERON_0000178',
        gene_id='ENSG00000000460',
        tmp_dir=tmp_dir,
        method='conditional'
        )

    return 0

if __name__ == '__main__':

    main()
