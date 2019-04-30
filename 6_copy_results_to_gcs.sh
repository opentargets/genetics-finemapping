#!/usr/bin/env bash
#

version_date=`date +%y%m%d`

# Copy results
gsutil -m rsync -r /home/ubuntu/results/finemapping/results gs://genetics-portal-staging/finemapping/$version_date

# Tar the logs and copy over
tar -zcvf logs.tar.gz /home/ubuntu/results/finemapping/logs
gsutil -m cp logs.tar.gz gs://genetics-portal-staging/finemapping/$version_date/logs.tar.gz
