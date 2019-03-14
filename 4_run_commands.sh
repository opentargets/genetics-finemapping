#!/usr/bin/env bash
#

set -euo pipefail

cores=4
instance_name="em-finemap-ssd"
instance_zone="europe-west1-d"

python 3_make_commands.py | parallel -j $cores

echo COMPLETE

# gcloud compute instances stop $instance_name --zone=$instance_zone