#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Splits the LD reference panel into multiple subfiles, to greatly speed up
# GCTA-cojo. Each subfile is a 3-Mb window, so that for any top_loci variant,
# we can use a subfile that has the top_loci variant +- 1 Mb.

import pandas as pd
import os
import subprocess as sp
import argparse


def main():

    # Parse args
    args = parse_args()

    window_size = int(3e6)
    window_spacing = int(1e6)

    chrom_lengths = pd.read_csv('/configs/grch38_chrom_lengths.tsv', sep='\t')

    # Loop through each chromosome and use plink to split the ref panel into
    # overlapping windows
    for index, row in chrom_lengths.iterrows():
        chr_ld_path = args.path.format(chrom=row['chrom'])
        print(chr_ld_path)

        window_start = int(0)
        while (window_start + window_size - window_spacing) < row['length']:
            # Define window and output path for subfile of main LD file
            window_end = window_start + window_size
            out_ld_path = chr_ld_path + '.{:d}_{:d}'.format(window_start, window_end)
            print("chr_ld_path:" + chr_ld_path)
            print("out_ld_path:" + out_ld_path)
            # plink requires a file to define the range to extract
            # We don't want a temp file, so we use bash process substitution <(...)
            range_str = ' '.join([str(row['chrom']), str(window_start), str(window_end), 'range1'])
            cmd = '/bin/bash 0_plink_extract.sh {} {} \'{}\''.format(chr_ld_path, out_ld_path, range_str)
            print(cmd)

            # Run plink
            os.system(cmd)
            cp = sp.run(cmd, shell=True, stderr=sp.STDOUT)
            if cp.returncode != 0:
                print('Failed on plink command:\n{}'.format(cmd))
                #return cp.returncode
            
            window_start = window_start + window_spacing

    return 0


def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--path',
                        metavar="<string>",
                        help='Path to LD reference; {chrom} in place of each chromosome name',
                        type=str,
                        required=True)
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    main()
