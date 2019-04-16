Fine-mapping pipeline
=====================

Todo:
- remove temp files option
- remove pos_37 pos_38


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
# Activate environment
source activate finemap

# Set spark paths
export PYSPARK_SUBMIT_ARGS="--driver-memory 8g pyspark-shell"
export SPARK_HOME=/home/ubuntu/software/spark-2.4.0-bin-hadoop2.7
export PYTHONPATH=$SPARK_HOME/python:$SPARK_HOME/python/lib/py4j-2.4.0-src.zip:$PYTHONPATH

# Create input manifest
python 1_make_input_file_manifest.py

# Execute workflow
bash execute_workflow.test.sh
```

Useful commands:

```
# Parse time taken for each run
grep "Time taken" logs/study_id=*/phenotype_id=*/bio_feature=*/chrom=*/logfile.txt
ls -rt logs/study_id=*/phenotype_id=*/bio_feature=*/chrom=*/logfile.txt | xargs grep "Time taken"

# List all
ls logs/study_id=*/phenotype_id=*/bio_feature=*/chrom=*/logfile.txt
```

# OLD

```
# Find output directories that don't contain the expected files
find ~/finemapping/output -mindepth 5 -maxdepth 5 -type d '!' -exec test -e "{}/top_loci.parquet" ';' -print
find ~/finemapping/output -mindepth 5 -maxdepth 5 -type d '!' -exec test -e "{}/credible_set.parquet" ';' -print

# List log files with errors
grep -n -i error ~/finemapping/logs/study_id\=*/cell_id\=*/group_id\=*/trait_id\=*/chrom\=*/logfile.txt | cut -f 1 -d ":" | sort | uniq

# How many have errors
grep -n -i error ~/finemapping/logs/study_id\=*/cell_id\=*/group_id\=*/trait_id\=*/chrom\=*/logfile.txt | cut -f 1 -d ":" | sort | uniq | wc -l

# Number of log files
ls ~/finemapping/logs/study_id\=*/cell_id\=*/group_id\=*/trait_id\=*/chrom\=*/logfile.txt | wc -l

# Count finished
find output -name "*.json" | wc -l
find output -name "logfile.txt" | wc -l

```

# TODO

- Only load required fields from parquet files
- Remove build 37 / 38 differentiation from scripts (still needed for excluding the MHC)
- Remove cell_id, trait_id, group_id.
- Keep study_id, phenotype_id, bio_feature
- Change loading of EAF, n_cases, n_total
- Currently fails for sex chromosomes
  - Need to replace X with 23 in plink file or when specifying gcta command
  - Need to impute sex in plink file for X for cojo to work
- Change outputs to json, not parquet
- Add type column to outputs
- Don't threshold the PP on credible set analysis

For readme:

- Manifest NAs must be represented with "None"
- Requirements for input files, e.g. if bio_feature not Null, it should be (Hive) paritioned by this field and chrom
- P-value threshold is specified in 1_scan_input_parquets.py. Set to 5e-8 for GWAS, and (0.05 / num_tests) for mol trait