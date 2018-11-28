Finemapping
===========

### Requirements
- GCTA (>= v1.91.3) must be available in `$PATH`
- [conda](https://conda.io/docs/)
- [Cromwell](https://cromwell.readthedocs.io/en/stable/) jar accessible via `$CROMWELL_JAR`
- git

### Setup enironment

```
git clone https://github.com/opentargets/finemapping.git
cd finemapping
bash setup.sh # Requires sudo log in part way through
conda env create -n finemapping --file environment.yaml
```

- Set `concurrent-job-limit` in `configs/cromwell.config`

### Usage

```
# Activate environment
source activate finemapping

# Start mysql docker
docker run -p 3306:3306 --name cromwell_myself -e MYSQL_ROOT_PASSWORD=cromwell_root_password -e MYSQL_DATABASE=cromwell_db -e MYSQL_USER=cromwell_user -e MYSQL_PASSWORD=cromwell_other_password -d mysql/mysql-server:5.5

# Create input manifest
python 1_make_input_file_manifest.py

# Execute workflow
bash execute_workflow.test.sh

```

# TODO
- Add logging
- Only load required fields from parquet files
- When loading, apply row filters in additions to row-group filters
- Add build 37 variant ID and position columns to sumstat files
- Add eaf_estimated column to the sumstat files and select this as eaf when loading the data (currently using MAF for molecular_qtl which is incorrect)
- Currently fails for sex chromosomes
  - Need to replace X with 23 in plink file or when specifying gcta command
  - Need to impute sex in plink file for X for cojo to work
- Only create output folders if there is anything to write
