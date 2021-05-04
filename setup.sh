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

# Install GCTA
cd $HOME
mkdir -p ~/software/gcta
cd ~/software/gcta
# Note that this URL may change - old versions aren't accessible at the same URL
wget https://cnsgenomics.com/software/gcta/bin/gcta_1.93.2beta.zip
unzip gcta_1.93.2beta.zip
cd gcta_1.93.2beta
echo export PATH="$PWD:\$PATH" >> ~/.profile
. ~/.profile

# Install plink
mkdir -p ~/software/plink
cd ~/software/plink
wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20201019.zip
unzip plink_linux_x86_64_20201019.zip
echo export PATH="$PWD:\$PATH" >> ~/.profile
. ~/.profile

# Install FINEMAP
mkdir -p ~/software/finemap
cd ~/software/finemap
wget http://www.christianbenner.com/finemap_v1.4_x86_64.tgz
tar -zxf finemap_v1.4_x86_64.tgz
ln -s finemap_v1.4_x86_64/finemap_v1.4_x86_64 finemap
sudo apt-get install libgomp1 # Not present by default it seems
echo export PATH="$PWD:\$PATH" >> ~/.profile
. ~/.profile

# Copy LD to VM
#gsutil -m cp ukb_v3_chr11* gs://genetics-portal-analysis/jeremy/ukb/tmp/
#gsutil -m cp gs://genetics-portal-analysis/jeremy/ukb/tmp/ukb_v3_chr11* ~/genetics-finemapping/data/ukb_downsampled10k/
mkdir -p ~/genetics-finemapping/data/ukb_downsampled10k/
gsutil -m cp -r gs://open-targets-ukbb/genotypes/ukb_v3_downsampled10k/ukb_v3_chr* ~/genetics-finemapping/data/ukb_downsampled10k/

# Install JRE
sudo apt install -yf openjdk-8-jre-headless openjdk-8-jdk
# sudo update-java-alternatives --list
# sudo update-java-alternatives --set java-1.8.0-openjdk-amd64

# Install parallel
sudo apt install -yf parallel

echo COMPLETE
