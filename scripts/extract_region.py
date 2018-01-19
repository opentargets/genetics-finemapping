#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#
"""
Extracts speicified region from GCTA formatted file

"""

import sys
import pandas as pd

def main():

    # Args
    insumstats = sys.argv[1]
    inbim = sys.argv[2]
    index_snp = sys.argv[3]
    window = int(sys.argv[4]) * 1000 # Take +-kb from SNP position
    min_maf = float(sys.argv[5])
    outf = sys.argv[6]

    # Load data
    data = pd.read_csv(insumstats, sep="\t")
    bim = pd.read_csv(inbim, sep="\t", header=None)
    bim.columns = ["chrom", "SNP", "x", "bp", "y", "z"]
    bim.index = bim.SNP

    # Keep only rows that are within speicified window around index_snp
    index_pos = bim.loc[index_snp, "bp"]
    in_window = ( (bim["bp"] >= (index_pos - window)) &
                  (bim["bp"] <= (index_pos + window)) )
    bim_subset = bim.loc[in_window, :]

    # Merge to sum stats
    merged = pd.merge(data, bim_subset, on="SNP", how="right")

    # Filter on maf
    to_keep = merged.freq.apply(freq_to_maf) >= min_maf
    merged = merged.loc[to_keep, :]

    # Make output df
    outdf = merged.rename(columns={"chrom":"Chr", "A1":"refA", "N":"n"})
    outdf = outdf.loc[:, ["Chr", "SNP", "bp", "refA", "freq", "b", "se", "p", "n"]]

    # Add required columns
    outdf["freq_geno"] = outdf["freq"]
    outdf["bC"] = outdf["b"]
    outdf["bC_se"] = outdf["se"]
    outdf["pC"] = outdf["p"]

    # Write
    outdf.to_csv(outf, sep="\t", index=None)


    return 0

def freq_to_maf(freq):
    freq = float(freq)
    if freq > 0.5:
        return 1 - freq
    else:
        return freq

if __name__ == '__main__':

    main()
