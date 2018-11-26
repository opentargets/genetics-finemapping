#!/usr/bin/env bash
#

set -euo pipefail

# Run cromwell
java -jar ~/software/cromwell_36/cromwell-36.jar run workflows/finemapping.wdl \
  --inputs configs/workflow.config.json

echo COMPLETE
