#!/usr/bin/env bash
#

#set -euo pipefail

cores="${CORES:-8}"
instance_name="em-finemapping-big"

python 3_make_commands.py | shuf | parallel -j $cores --bar --joblog parallel.jobs.log

echo COMPLETE
