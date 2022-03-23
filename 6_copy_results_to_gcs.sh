#!/usr/bin/env bash
#

version_date=$1

# Copy results
gsutil -m rsync -r $HOME/genetics-finemapping/results/ gs://genetics-portal-dev-staging/finemapping/$version_date

# Tar the logs and copy over
# This can take a very long time, so you may not want to keep the logs at all
#tar -zcvf logs.tar.gz $HOME/genetics-finemapping/logs
#gsutil -m cp logs.tar.gz gs://genetics-portal-dev-staging/finemapping/$version_date/logs.tar.gz
