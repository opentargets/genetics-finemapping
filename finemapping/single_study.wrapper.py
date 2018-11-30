#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import sys
import os
import argparse
import logging
import main as fm
import pprint
from datetime import datetime
import yaml

def main():

    # Args
    args = parse_args()
    start_time = datetime.now()

    # Make output and temp directories
    os.makedirs(os.path.split(args.toploci)[0], exist_ok=True)
    os.makedirs(os.path.split(args.credset)[0], exist_ok=True)
    os.makedirs(os.path.split(args.log)[0], exist_ok=True)
    os.makedirs(args.tmpdir , exist_ok=True)

    # Start logging
    logger = make_logger(args.log)
    logger.info('Started finemapping pipeline')
    logger.info('Printing args: \n' + pprint.pformat(vars(args), indent=2))

    # Load analysis config file
    with open(args.config_file, 'r') as in_h:
        config_dict = yaml.load(in_h)
    logger.info('Analysis config: \n' + pprint.pformat(config_dict, indent=2))

    # Run
    top_loci, credset_results = fm.run_single_study(
        in_pq=args.pq,
        in_plink=args.ld,
        study_id=args.study_id,
        cell_id=args.cell_id,
        group_id=args.group_id,
        trait_id=args.trait_id,
        chrom=args.chrom,
        analysis_config=config_dict,
        tmp_dir=args.tmpdir,
        method=args.method,
        logger=logger
    )

    # Write output
    logger.info('Writing outputs')
    # if top_loci.shape[0] > 0:
    top_loci.to_parquet(
        fname=args.toploci,
        engine='fastparquet',
        compression='snappy'
    )
    # if credset_results is not None:
    credset_results.to_parquet(
        fname=args.credset,
        engine='fastparquet',
        compression='snappy'
    )

    # Log time taken
    logger.info('Time taken: {0}'.format(
        datetime.now() - start_time
    ))
    logger.info('Finished!')

    return 0

def parse_args():
    ''' Load command line args
    '''
    p = argparse.ArgumentParser()

    # Add input files
    p.add_argument('--pq',
                   metavar="<file>",
                   help=("Input: parquet file containing summary stats"),
                   type=str,
                   required=True)
    p.add_argument('--ld',
                   metavar="<file>",
                   help=("Input: plink file to estimate LD from"),
                   type=str,
                   required=True)
    p.add_argument('--config_file',
                   metavar="<file>",
                   help=("Input: analysis config file"),
                   type=str,
                   required=True)

    # Add study identifier args
    for key in ['study_id', 'trait_id', 'chrom']:
        p.add_argument('--{0}'.format(key),
                       metavar="<str>",
                       help=("{0} to extract from pq".format(key)),
                       type=str,
                       required=True)
    # Add mol trait study identifier args
    for key in ['cell_id', 'group_id']:
        p.add_argument('--{0}'.format(key),
                       metavar="<str>",
                       help=("{0} to extract from pq".format(key)),
                       type=str,
                       required=False)

    # Add methological args
    p.add_argument('--method',
                   metavar="[conditional|distance]",
                   help=("Which method to run, either with conditional analysis"
                         " (gcta-cojo) or distance based with conditional"),
                   type=str,
                   choices=['conditional', 'distance'],
                   required=True)

    # Add output file
    p.add_argument('--toploci',
                   metavar="<file>",
                   help=("Output: top loci parquet file"),
                   type=str,
                   required=True)
    p.add_argument('--credset',
                   metavar="<file>",
                   help=("Output: credible set parquet file"),
                   type=str,
                   required=True)
    p.add_argument('--log',
                   metavar="<file>",
                   help=("Output: log file"),
                   type=str,
                   required=True)
    p.add_argument('--tmpdir',
                   metavar="<file>",
                   help=("Output: temp dir"),
                   type=str,
                   required=True)


    args = p.parse_args()

    # Convert "None" strings to None type
    for arg in vars(args):
        if getattr(args, arg) == "None":
            setattr(args, arg, None)

    return args

def make_logger(log_file):
    ''' Creates a logging handle.
    '''
    # Basic setup
    logging.basicConfig(
        level=logging.DEBUG,
        datefmt='%Y-%m-%d %H:%M:%S',
        stream=None)
    # Create formatter
    logFormatter = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
    rootLogger = logging.getLogger(__name__)
    # Add file logging
    fileHandler = logging.FileHandler(log_file, mode='w')
    fileHandler.setFormatter(logFormatter)
    rootLogger.addHandler(fileHandler)
    # Add stdout logging
    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    rootLogger.addHandler(consoleHandler)
     # Prevent logging from propagating to the root logger
    rootLogger.propagate = 0

    return rootLogger

if __name__ == '__main__':

    main()
