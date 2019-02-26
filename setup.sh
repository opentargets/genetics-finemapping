#!/usr/bin/env bash
#
# Must be ran as root

# set -euo pipefail

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

# Install docker: https://docs.docker.com/install/linux/docker-ce/ubuntu/#install-using-the-repository
cd $HOME
sudo apt-get remove docker docker-engine docker.io
sudo apt-get update
sudo apt-get install \
    apt-transport-https \
    ca-certificates \
    curl \
    software-properties-common
sudo curl -fsSL https://download.docker.com/linux/ubuntu/gpg | apt-key add -
sudo apt-key fingerprint 0EBFCD88
sudo add-apt-repository \
   "deb [arch=amd64] https://download.docker.com/linux/ubuntu \
   $(lsb_release -cs) \
   stable"
sudo apt-get update
sudo apt-get install docker-ce
sudo apt install docker-compose
# Create docker group and add USER
sudo groupadd docker
sudo usermod -aG docker $USER
# Run test
docker run hello-world # Test install

# Install JRE
sudo apt install -yf openjdk-8-jre-headless openjdk-8-jdk
# sudo update-java-alternatives --list
# sudo update-java-alternatives --set java-1.8.0-openjdk-amd64


echo COMPLETE
