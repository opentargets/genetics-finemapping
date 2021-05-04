#!/usr/bin/env bash
#

version_date=`date +%y%m%d`

# Copy results
gsutil -m rsync -r $HOME/genetics-finemapping/results/ gs://genetics-portal-dev-staging/finemapping/$version_date

# Tar the logs and copy over
tar -zcvf logs.tar.gz $HOME/genetics-finemapping/logs
gsutil -m cp logs.tar.gz gs://genetics-portal-dev-staging/finemapping/$version_date/logs.tar.gz
