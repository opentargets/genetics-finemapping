#!/usr/bin/env bash
#

set -euo pipefail

cores=59
instance_name="em-finemapping-big"

python 3_make_commands.py | shuf | parallel -j $cores

echo COMPLETE

openstack server suspend $instance_name
