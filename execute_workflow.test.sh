#!/usr/bin/env bash
#

set -euo pipefail

CROMWELL_JAR=/home/emountjoy_statgen/software/cromwell/cromwell-38.jar

# Run cromwell
mkdir -p logs
java -Dconfig.file=configs/cromwell.config \
     -jar $CROMWELL_JAR run workflows/finemapping.wdl \
     --inputs configs/workflow.config.json \
     > cromwell_log.txt

echo COMPLETE
