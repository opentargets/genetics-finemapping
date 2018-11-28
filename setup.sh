#!/usr/bin/env bash
#
# Must be ran as root

set -euo pipefail

# Install cromwell
cd $HOME
mkdir -p software/cromwell
cd software/cromwell
wget https://github.com/broadinstitute/cromwell/releases/download/36/cromwell-36.jar
wget https://github.com/broadinstitute/cromwell/releases/download/36/womtool-36.jar
echo export export CROMWELL_JAR=$PWD/cromwell-36.jar >> ~/.profile
. ~/.profile

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
wget https://cnsgenomics.com/software/gcta/gcta_1.91.7beta.zip
unzip gcta_1.91.7beta.zip
cd gcta_1.91.7beta
echo export PATH="$PWD:\$PATH" >> ~/.profile
. ~/.profile

# Log in as superuser
cd $HOME
sudo su -

# Install docker: https://docs.docker.com/install/linux/docker-ce/ubuntu/#install-using-the-repository
# apt-get remove docker docker-engine docker.io
# apt-get update
# apt-get install \
#     apt-transport-https \
#     ca-certificates \
#     curl \
#     software-properties-common
# curl -fsSL https://download.docker.com/linux/ubuntu/gpg | apt-key add -
# apt-key fingerprint 0EBFCD88
# add-apt-repository \
#    "deb [arch=amd64] https://download.docker.com/linux/ubuntu \
#    $(lsb_release -cs) \
#    stable"
# apt-get update
# apt-get install docker-ce
# docker run hello-world # Test install

# Install JRE
apt install -yf openjdk-8-jre-headless

echo COMPLETE
