#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import json

def main():

    data = [
        {
            'type': 'gwas',
            'study_id': 'NEALEUKB_50',
            'cell_id': None,
            'group_id': None,
            'gene_id': None,
            'trait_id': 'UKB_50',
            'chrom': '22',
            'method': 'conditional',
            'ld_ref': 'input/ld/EUR.{chrom}.1000Gp3.20130502'
        },
        {
            'type': 'gwas',
            'study_id': 'NEALEUKB_50',
            'cell_id': None,
            'group_id': None,
            'gene_id': None,
            'trait_id': 'UKB_50',
            'chrom': '21',
            'method': 'conditional',
            'ld_ref': 'input/ld/EUR.{chrom}.1000Gp3.20130502'
        },
        {
            'type': 'molecular_qtl',
            'study_id': 'GTEX7',
            'cell_id': 'UBERON_0000178',
            'group_id': 'ENSG00000000460',
            'gene_id': 'ENSG00000000460',
            'trait_id': 'eqtl',
            'chrom': '1',
            'method': 'conditional',
            'ld_ref': 'input/ld/EUR.{chrom}.1000Gp3.20130502'
        }
    ]

    with open('manifest.json', 'w') as out_h:
        json.dump(data, out_h, indent=2)

    return 0

if __name__ == '__main__':

    main()
