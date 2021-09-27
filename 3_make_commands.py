#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#
# Reads the manifest file and makes commands
#

import os
import sys
import json
import argparse
import gzip
import yaml

def main():

    # Args
    args = parse_args()
    in_manifest = 'configs/manifest.json.gz'
    out_todo = 'commands_todo.txt.gz'
    out_done = 'commands_done.txt.gz'

    # Pipeline args
    script = 'finemapping/single_study.wrapper.py'
    analysis_config = 'configs/analysis.config.yaml'

    with open(analysis_config, 'r') as in_h:
        config_dict = yaml.load(in_h, Loader=yaml.FullLoader)

    # Open command files
    todo_h = gzip.open(out_todo, 'w')
    done_h = gzip.open(out_done, 'w')
    
    # Iterate over manifest
    with gzip.open(in_manifest, 'r') as in_mani:
        for line in in_mani:

            # Parse
            rec = json.loads(line.decode().rstrip())

            # Build command
            cmd = [
                'python',
                os.path.abspath(script),
                '--pq', os.path.abspath(rec['in_pq']),
                '--ld', os.path.abspath(rec['in_ld']),
                '--config_file', os.path.abspath(analysis_config),
                '--study_id', rec['study_id'],
                '--phenotype_id', rec['phenotype_id'],
                '--bio_feature', rec['bio_feature'],
                '--type', rec['type'],
                '--chrom', rec['chrom'],
                '--method', rec['method'],
                '--pval_threshold', rec['pval_threshold'],
                '--toploci', os.path.abspath(rec['out_top_loci']),
                '--credset', os.path.abspath(rec['out_credset']),
                '--finemap', os.path.abspath(rec['out_finemap']),
                '--tmpdir', os.path.abspath(rec['tmpdir']),
                '--log', os.path.abspath(rec['out_log']),
                '--delete_tmpdir'
            ]
            cmd_str = ' '.join([str(arg) for arg in cmd])

            # Skip if both toploci and credset outputs exist
            if (os.path.exists(rec['out_top_loci']) and
                    os.path.exists(rec['out_credset']) and
                    (not config_dict['run_finemap'] or os.path.exists(rec['out_finemap']))):
                done_h.write((cmd_str + '\n').encode())
                continue
            else:
                todo_h.write((cmd_str + '\n').encode())
                if not args.quiet:
                    print(cmd_str)
    
    # Close files
    done_h.close()
    todo_h.close()

    return 0

def parse_args():
    ''' Load command line args
    '''
    p = argparse.ArgumentParser()

    # Add input files
    p.add_argument('--quiet',
                   help=("Don't print commands to stdout"),
                   action='store_true')

    args = p.parse_args()
    return args

if __name__ == '__main__':

    main()
