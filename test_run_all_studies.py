#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import finemapping.ot_pipeline
import os
import yaml

def main():

    # Load pipleine config file
    config_file = 'configs/pipeline.config.yaml'
    with open(config_file, 'r') as in_h:
        pipe_config = yaml.load(in_h)

    print('Pipe config:', pipe_config)

    # Run
    finemapping.ot_pipeline.run_all_studies(**pipe_config)

    return 0

if __name__ == '__main__':

    main()
