#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import finemapping.ot_pipeline
import os

def main():


    # Set up test params
    tmp_dir = 'tmp/NEALEUKB_50'
    out_dir = 'output/NEALEUKB_50'

    # Create output tmp dir
    # os.makedirs(tmp_dir, exist_ok=True)

    # Run example
    finemapping.ot_pipeline.run_single_study(
        in_pq='input/gwas/NEALEUKB_50',
        in_plink='input/ld/EUR.{chrom}.1000Gp3.20130502',
        study_id='NEALEUKB_50',
        config_file='config/config.yaml',
        chrom='22',
        tmp_dir=tmp_dir,
        # method='distance'
        method='conditional'
        )


    return 0

if __name__ == '__main__':

    main()
