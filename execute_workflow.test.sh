#!/usr/bin/env bash
#

set -euo pipefail

CROMWELL_JAR=/Users/em21/software/cromwell/cromwell-36.jar

# Run cromwell
mkdir -p logs
java -Dconfig.file=configs/cromwell.local.config \
     -jar $CROMWELL_JAR run workflows/finemapping.wdl \
     --inputs configs/workflow.config.json \
     > cromwell_log.txt

echo COMPLETE
