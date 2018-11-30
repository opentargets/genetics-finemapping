Finemapping
===========

### Requirements
- GCTA (>= v1.91.3) must be available in `$PATH`
- [conda](https://conda.io/docs/)
- [Cromwell](https://cromwell.readthedocs.io/en/stable/) jar accessible via `$CROMWELL_JAR`
- git

### Setup environment

```
git clone https://github.com/opentargets/finemapping.git
cd finemapping
bash setup.sh
conda env create -n finemapping --file environment.yaml
```

### Edit config files

- `analysis.config.yaml`: Fine-mapping analysis parameters
- `configs/input_files.config.tsv`: Manifest file specifying all studies to run fine-mapping on. File should be tab-separated with no header and the following column order:
  * `in_pq`: Input parquet file
  * `ld_ref`: Pattern of plink ld ref file. Use "{chrom}" as wildcard.
  * `study_id`: Study ID in `in_pq`
  * `cell_id`: Cell ID in `in_pq`
  * `group_id`: Group ID in `in_pq`
  * `trait_id`: Trait ID in `in_pq`
  * `chrom`: Chrom in `in_pq`
  * `method`: Method either "conditional" or "distance" based
  * `out_top_loci`: Output top loci parquet
  * `out_credset`: Output credible set parquet
  * `out_log`: Output log file
  * `tmpdir`: Temp directory
- `configs/workflow.config.json`: Input files for the Cromwell workflow
- `configs/cromwell.config`: Cromwell configuration file
  * recommended: set `concurrent-job-limit`

### Usage

```
# Start mysql server docker
mkdir -p mysql/data
mkdir -p mysql/init
docker run \
  -v /home/ubuntu/finemapping/mysql/data:/var/lib/mysql \
  -v /home/ubuntu/finemapping/mysql/init:/docker-entrypoint-initdb.d \
  -p 3306:3306 \
  --name cromwell_myself \
  -e MYSQL_ROOT_PASSWORD=cromwell_root_password \
  -e MYSQL_DATABASE=cromwell_db \
  -e MYSQL_USER=cromwell_user \
  -e MYSQL_PASSWORD=cromwell_other_password \
  -d mysql/mysql-server:5.5

# Activate environment
source activate finemapping

# Create input manifest
python 1_make_input_file_manifest.py

# Execute workflow
bash execute_workflow.test.sh
```

# TODO
- Only load required fields from parquet files
- Add build 37 variant ID and position columns to sumstat files
- Add eaf_estimated column to the sumstat files and select this as eaf when loading the data (currently using MAF for molecular_qtl which is incorrect)
- Currently fails for sex chromosomes
  - Need to replace X with 23 in plink file or when specifying gcta command
  - Need to impute sex in plink file for X for cojo to work
- Only create output folders if there is anything to write
- Test if no top loci are found or not credible sets are found
