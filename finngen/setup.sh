#!/usr/bin/env bash
#
# Must be run as root

# set -euo pipefail

sudo apt-get update
sudo apt-get install unzip

# Install conda
cd $HOME
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b -p $HOME/miniconda
echo export PATH="$HOME/miniconda/bin:\$PATH" >> ~/.profile
. ~/.profile

sudo apt install -yf parallel
