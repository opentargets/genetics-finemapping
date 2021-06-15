#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Jeremy Schwartzentruber
#
# Reads the manifest file and makes commands
#

import os
import sys
import json
import argparse
import gzip

def main():

    # Args
    args = parse_args()
    in_manifest = 'finngen.manifest.json'
    out_todo = 'commands_todo.txt.gz'
    out_done = 'commands_done.txt.gz'

    # Pipeline args
    script = 'finngen_finemapping_ingest.py'
    analysis_config = '../configs/analysis.config.yaml'

    # Open command files
    todo_h = gzip.open(out_todo, 'w')
    done_h = gzip.open(out_done, 'w')
    
    # Iterate over manifest
    with open(in_manifest, 'r') as in_mani:
        for line in in_mani:

            # Parse
            rec = json.loads(line.rstrip())

            # Build command
            cmd = [
                'python',
                os.path.abspath(script),
                '--config_file', os.path.abspath(analysis_config),
                '--study_id', rec['study_id'],
                '--in_snp', os.path.abspath(rec['in_snp']),
                '--toploci', os.path.abspath(rec['out_top_loci']),
                '--credset', os.path.abspath(rec['out_credset']),
                '--log', os.path.abspath(rec['out_log']),
            ]
            cmd_str = ' '.join([str(arg) for arg in cmd])

            # Skip if both toploci and credset outputs exist
            if (os.path.exists(rec['out_top_loci']) and
                    os.path.exists(rec['out_credset'])):
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
