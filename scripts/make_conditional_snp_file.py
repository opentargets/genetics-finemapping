#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#
"""
Takes the output of multi loci clustering script and outputs a file containing
snps on which to condition.
"""

import sys

def main():

    # Args
    inf = sys.argv[1]
    outf = sys.argv[2]
    index_snp = sys.argv[3]

    # Find conditional snps
    conditional_snps = None
    with open(inf, "r") as in_h:
        in_h.readline()
        for line in in_h:
            parts = line.rstrip().split()
            if parts[0] == index_snp:
                try:
                    conditional_snps = parts[1].split(";")
                except IndexError:
                    conditional_snps = []

    # Raise error if none are found
    if conditional_snps is None:
        sys.exit("Error: {0} not found in {1}".format(index_snp, inf))

    # Write output
    with open(outf, "w") as out_h:
        for conditional_snp in conditional_snps:
            out_h.write(conditional_snp + "\n")

    return 0

if __name__ == '__main__':

    main()
