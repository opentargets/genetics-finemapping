#!/usr/bin/env bash
#
# Must be ran as root

# set -euo pipefail

sudo apt-get update
sudo apt-get install unzip

# Install conda
cd $HOME
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b -p $HOME/miniconda
echo export PATH="$HOME/miniconda/bin:\$PATH" >> ~/.profile
. ~/.profile

# Install GCTA
cd $HOME
mkdir -p software/gcta
cd software/gcta
wget https://cnsgenomics.com/software/gcta/gcta_1.92.0beta3.zip
unzip gcta_1.92.0beta3.zip
cd gcta_1.92.0beta3
echo export PATH="$PWD:\$PATH" >> ~/.profile
. ~/.profile

# Install JRE
sudo apt install -yf openjdk-8-jre-headless openjdk-8-jdk
# sudo update-java-alternatives --list
# sudo update-java-alternatives --set java-1.8.0-openjdk-amd64

# Install parallel
sudo apt install -yf parallel

echo COMPLETE
