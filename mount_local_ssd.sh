#!/usr/bin/env bash
#

set -euo pipefail

lsblk
sudo mkfs.ext4 -F /dev/sdb
sudo mkdir -p ~/scratch
sudo mount /dev/sdb ~/scratch
sudo chmod a+w ~/scratch/

echo COMPLETE
