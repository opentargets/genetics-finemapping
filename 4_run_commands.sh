#!/usr/bin/env bash
#

set -euo pipefail

cores=8
instance_name="em-finemapping-big"

python 3_make_commands.py | shuf | parallel -j $cores --bar --joblog logs/parallel.jobs.log

echo COMPLETE

#openstack server suspend $instance_name
