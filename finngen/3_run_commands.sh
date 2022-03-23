#!/usr/bin/env bash
#
NCORES=$1

set -euo pipefail

mkdir -p logs

python 2_make_commands.py | shuf | parallel -j $NCORES --bar --joblog logs/parallel.jobs.log

echo COMPLETE

#openstack server suspend $instance_name
