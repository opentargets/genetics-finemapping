#!/usr/bin/env bash
#

set -euo pipefail

cores=7
mkdir -p logs

python 2_make_commands.py | shuf | parallel -j $cores --bar --joblog logs/parallel.jobs.log

echo COMPLETE

#openstack server suspend $instance_name
