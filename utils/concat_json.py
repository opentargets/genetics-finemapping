#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import argparse
import gzip

def main():

    # Args
    args = parse_args()

    # Concat together
    with gzip.open(args.out, 'w') as out_h:
        for inf in args.in_json:
            with gzip.open(inf, 'r') as in_h:
                for line in in_h:
                    line = line.decode().rstrip()
                    out_h.write((line + '\n').encode())

    return 0

def parse_args():
    ''' Load command line args
    '''
    p = argparse.ArgumentParser()

    # Add input files
    p.add_argument('--in_json',
                   metavar="<file>",
                   help=("List of json files to concatenate"),
                   type=str,
                   nargs='+',
                   required=True)
    p.add_argument('--out',
                   metavar="<file>",
                   help=("Concatenated json file"),
                   type=str,
                   required=True)

    args = p.parse_args()

    return args

if __name__ == '__main__':

    main()

