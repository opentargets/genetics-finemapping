#!/usr/bin/env bash
#

set -euo pipefail

cores=7
python 3_make_commands.py | shuf | parallel -j $cores --bar --joblog logs/parallel.jobs.log

echo COMPLETE
