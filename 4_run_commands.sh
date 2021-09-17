#!/usr/bin/env bash
#
NCORES=$1

set -euo pipefail

python 3_make_commands.py | shuf | parallel -j $NCORES --bar --joblog logs/parallel.jobs.log

echo COMPLETE
