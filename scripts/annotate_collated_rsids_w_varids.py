#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import sys
import os
import pandas as pd
import argparse

def main():

    # Args
    args = parse_args()

    # Load cred sets
    cred = pd.read_csv(args.inf, sep="\t", header=0)

    # Create set of all rsids
    rsid_set = set(cred["locus_index_snp"]) | set(cred["SNP"])

    # Load lines from sumstat file
    iter_csv = pd.read_csv(args.sumstat, sep="\t", header=0, iterator=True, chunksize=100000)
    sumstat = pd.concat([chunk[chunk['rsid'].isin(rsid_set), ["rsid", "snpid"]] for chunk in iter_csv])

    # Make map dict
    rsid_varid = dict(zip(sumstat["rsid"], sumstat["snpid"]))
    locus_index_varid = cred["locus_index_snp"].apply(lambda rsid: rsid_varid[rsid])
    varid = cred["SNP"].apply(lambda rsid: rsid_varid[rsid])

    # Insert var ids
    cred.insert(2, "locus_index_varid", locus_index_varid)
    cred.insert(5, "varid", varid)

    # Save
    cred.to_csv(args.outf, sep="\t", index=None, compression="gzip")

    return 0

def parse_args():
    """ Load command line args.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--inf', metavar="<file>", help=('Input credible set file'), type=str, required=False)
    parser.add_argument('--outf', metavar="<file>", help=('Output file'), type=str, required=False)
    parser.add_argument('--sumstat', metavar="<file>", help=('Input reference sumstat file from which to load rsid to variant ID mappings'), type=str, required=False)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
