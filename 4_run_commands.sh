#!/usr/bin/env bash
#

set -euo pipefail

cores=55
instance_name="em-finemapping-big"

python 3_make_commands.py | parallel -j $cores

echo COMPLETE

# openstack server suspend $instance_name
