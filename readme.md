Open Target Genetics fine-mapping pipeline
==========================================

Fine-mapping pipeline for Open Targets Genetics. In brief, the method is:
1. Detect independent loci across the summary stat file using either (i) GCTA-cojo and a given plink file as an LD reference, (ii) distance based clumping. Method specified with `--method` argument.
2. If `--method conditional`, for each independent locus condition on all other surrounding loci (configurable with `cojo_wind`).
3. Perform approximate Bayes factor credible set analysis for each independent locus.

For FinnGen, we incorporate each new release by directly taking the SuSIE fine-mapping outputs from FinnGen to determine top loci. Each time this is done, we have to be careful to NOT to run GCTA fine-mapping on FinnGen sumstats, since their SuSIE fine-mapping is superior (and we don't have a good LD reference).

### Requirements
- Spark v2.4.0
- GCTA (>= v1.91.3) must be available in `$PATH`
- [conda](https://conda.io/docs/)
- GNU parallel

### Setup environment

```
git clone https://github.com/opentargets/genetics-finemapping.git
cd ~/genetics-finemapping
bash setup.sh
. ~/.profile # Reload profile so that conda works
conda env create -n finemap --file environment.yaml
```

### Configure pipeline

Many of the pipeline parameters must first be specified in the analysis config file: `configs/analysis.config.yaml`

### Run a single study

A single study can be fine-mapped using the single study wrapper

```
# Activate environment
source activate finemap

# Edit config file (this needs selecting with --config_file arg)
nano configs/analysis.config.yaml

# View args
$ python finemapping/single_study.wrapper.py --help
usage: single_study.wrapper.py [-h] --pq <file> --ld <file> --config_file
                               <file> --type <str> --study_id <str> --chrom
                               <str> [--phenotype_id <str>]
                               [--bio_feature <str>] --method
                               [conditional|distance] --pval_threshold <float>
                               --toploci <file> --credset <file> --log <file>
                               --tmpdir <file> [--delete_tmpdir]

optional arguments:
  -h, --help            show this help message and exit
  --pq <file>           Input: parquet file containing summary stats
  --ld <file>           Input: plink file to estimate LD from
  --config_file <file>  Input: analysis config file
  --type <str>          type to extract from pq
  --study_id <str>      study_id to extract from pq
  --chrom <str>         chrom to extract from pq
  --phenotype_id <str>  phenotype_id to extract from pq
  --bio_feature <str>   bio_feature to extract from pq
  --method [conditional|distance]
                        Which method to run, either with conditional analysis
                        (gcta-cojo) or distance based with conditional
  --pval_threshold <float>
                        P-value threshold to be considered "significant"
  --run_finemap         If True, then run FINEMAP
  --toploci <file>      Output: top loci json file
  --credset <file>      Output: credible set json file
  --finemap <file>      Output: finemap snp probabilities file
  --log <file>          Output: log file
  --tmpdir <file>       Output: temp dir
  --delete_tmpdir       Remove temp dir when complete

Note: The capability of running FINEMAP has been used but not extensively tested.
```

### Running the pipeline

#### Step 1: Prepare input data

[Prepare summary statistic files](https://github.com/opentargets/genetics-sumstat-data/), including ["significant" window extraction](https://github.com/opentargets/genetics-sumstat-data/tree/master/filters/significant_window_extraction) to reduce input size.

Prepare LD references in plink `bed|bim|fam` format, currently using [UK Biobank downsampled to 10K individuals and lifted over to GRCh38](https://github.com/opentargets/genetics-backend/tree/master/reference_data/uk_biobank_v3).

Download to local machine.
To process only new data, we should only download data from new "significant window" directories.
```
cd ~/genetics-finemapping
mkdir -p $HOME/genetics-finemapping/data/filtered/significant_window_2mb/gwas/
mkdir -p $HOME/genetics-finemapping/data/filtered/significant_window_2mb/molecular_trait/

gsutil -m rsync -r gs://genetics-portal-dev-sumstats/filtered/significant_window_2mb/gwas/ $HOME/genetics-finemapping/data/filtered/significant_window_2mb/gwas/
gsutil -m rsync -r gs://genetics-portal-dev-sumstats/filtered/significant_window_2mb/molecular_trait/ $HOME/genetics-finemapping/data/filtered/significant_window_2mb/molecular_trait/

# Remove FinnGen GWAS, since we don't run fine-mapping for them
rm -r $HOME/genetics-finemapping/data/filtered/significant_window_2mb/gwas/FINNGEN*

# Remove extra files that can screw up data loading later
find /home/js29/genetics-finemapping/data/filtered/significant_window_2mb -name "*_SUCCESS" | wc -l
find /home/js29/genetics-finemapping/data/filtered/significant_window_2mb -name "*_SUCCESS" -delete
```

#### Step 2: Prepare environment

```
# Set spark paths
export PYSPARK_SUBMIT_ARGS="--driver-memory 80g pyspark-shell"
export SPARK_HOME=/home/ubuntu/software/spark-2.4.0-bin-hadoop2.7
export PYTHONPATH=$SPARK_HOME/python:$SPARK_HOME/python/lib/py4j-2.4.0-src.zip:$PYTHONPATH
```

#### Step 3: Make manifest file

The manifest file specifies all analyses to be run. The manifest is a JSON lines file with each line containing the following fields:

```json
{
  "type": "gwas",
  "study_id": "NEALE2_50_raw",
  "phenotype_id": null,
  "bio_feature": null,
  "chrom": "6",
  "in_pq": "/home/ubuntu/data/sumstats/filtered/significant_window_2mb/gwas/NEALE2_50_raw.parquet",
  "in_ld": "/home/ubuntu/data/genotypes/ukb_v3_downsampled10k_plink/ukb_v3_chr{chrom}.downsampled10k",
  "out_top_loci": "/home/ubuntu/results/finemapping/output/study_id=NEALE2_50_raw/phenotype_id=None/bio_feature=None/chrom=6/top_loci.json.gz",
  "out_credset": "/home/ubuntu/results/finemapping/output/study_id=NEALE2_50_raw/phenotype_id=None/bio_feature=None/chrom=6/credible_set.json.gz",
  "out_log": "/home/ubuntu/results/finemapping/logs/study_id=NEALE2_50_raw/phenotype_id=None/bio_feature=None/chrom=6/logfile.txt",
  "tmpdir": "/home/ubuntu/results/finemapping/tmp/study_id=NEALE2_50_raw/phenotype_id=None/bio_feature=None/chrom=6",
  "method": "conditional",
  "pval_threshold": 5e-08
}
```

Note that the wildcard `{chrom}` can be used in `in_ld` field.

The manifest file can be automatically generated using:

```
cd ~/genetics-finemapping

# Edit the Args and Paths in `1_scan_input_parquets.py`
nano 1_scan_input_parquets.py

# Reads variants filtered for p value, and outputs a single json record in
# tmp/filtered_input for each study/chromosome combination with at least one
# significant variant. Takes a couple of minutes for 200 GWAS.
time python 1_scan_input_parquets.py

# Creates manifest file, one job per study/chromosome. Output path `configs/manifest.json.gz`
python 2_make_manifest.py
```

#### Step 4: Run pipeline

```
mkdir logs
tmux   # So run continues if connection is lost

# Edit args in `4_run_commands.sh` (e.g. number of cores) and then
NCORES=30
time bash 4_run_commands.sh $NCORES
zcat commands_todo.txt.gz | shuf | parallel -j $NCORES --bar --joblog logs/parallel.jobs2.log

# Exit tmux with Ctrl+b then d
```

The above command will run all analyses specified in the manifest using GNU parallel. It will create two files `commands_todo.txt.gz` and `commands_done.txt.gz` showing which analyses have not yet/already been done. The pipeline can be stopped at any time and restarted without repeating any completed analyses. You can safely regenerate the `commands_*.txt.gz` commands whilst the pipeline is running using `python 3_make_commands.py --quiet`.

If you get this error:
  ModuleNotFoundError: No module named 'dask'
then I've solved it just by deactivating conda and reactivating. This seems to happen especially when using tmux... I'm not sure why.

#### Step 5: Process the results

```
rm -r /home/js29/genetics-finemapping/tmp/*

# Combine the results of all the individual analyses
# This step can be slow/inefficient due to Hadoop many small files problem
# You are likely to get out of memory errors if you don't increase the java
# heap space available to Spark with PYSPARK_SUBMIT_ARGS.
export PYSPARK_SUBMIT_ARGS="--driver-memory 80g pyspark-shell"
time python 5_combine_results.py
```

```
# The below steps were used when we found duplicate top_loci, which was due
# to duplicated lines in the eQTL catalogue ingest. This has since been fixed,
# so the below should not be needed.
# Concatenate together all top_loci and credset files
time find output -name "top_loci.json.gz" | while read -r file; do zcat -f "$file"; done | gzip > top_loci.concat.json.gz &
time find output -name "credible_set.json.gz" | while read -r file; do zcat -f "$file"; done | gzip > credible_set.concat.json.gz

# Remove duplicates
# This should only be necessary because when we last ingested eQTL catalogue
# I failed to remove duplicate rows first.
time zcat top_loci.concat.json.gz | sort | uniq | gzip > top_loci.dedup.json.gz &
time zcat credible_set.concat.json.gz | sort | uniq | gzip > credible_set.dedup.json.gz

time python 5_combine_results_rmdup.py
```

```
# Make a note as to what this finemapping run contained. E.g.:
echo "Run with updated QTL datasets, and updated GWAS catalog studies. Re-ran all previous studies, since QTL datasets are the bulk of the fine-mapping work. Fixed an issue with 210825 version." > results/README.txt

# Copy the results to GCS
bash 6_copy_results_to_gcs.sh
```

Number of top_loci raw: 1,623,534
Number of top_loci after dups removed: 1,541,938
Number of top_loci in previous version: 689,726

Number of credset rows raw: 44,180,640
Number of credset rows after dups removed: 40,910,064


#### Step 6: Merge with previous fine-mapping results

Steps like the below are needed if we are adding on to previous fine-mapping results, rather than recomputing everything. This is also needed each time to incorporate FinnGen.
```
mkdir -p finemapping_temp/210923
mkdir -p finemapping_to_merge/finngen_210515
mkdir -p finemapping_merged
gsutil -m rsync -r gs://genetics-portal-dev-staging/finemapping/210923 finemapping_temp/210923
gsutil -m rsync -r gs://genetics-portal-dev-staging/finemapping/finngen_210515 finemapping_to_merge/finngen_210515

# Filter out any FinnGen loci from last release
zcat finemapping_temp/210923/top_loci.json.gz | grep -v "FINNGEN" | gzip > finemapping_temp/210923/top_loci_no_FINNGEN.json.gz

zcat finemapping_temp/210923/top_loci_no_FINNGEN.json.gz \
    finemapping_to_merge/finngen_210515/top_loci.json.gz \
    | gzip > finemapping_merged/top_loci.json.gz

python 7a_credset_remove_finngen.py # Write new dataset without FinnGen

zcat finemapping_temp/210923/credset/part*.json.gz | wc -l
zcat finemapping_to_merge/210923/credset/part*.json.gz | wc -l # After FinnGen loci removed

# Merge all non-FinnGen credsets with latest FinnGen credsets
python 7_merge_finemap_results.py

zcat finemapping_merged/credset/part*.json.gz | wc -l

gsutil -m rsync -r $HOME/genetics-finemapping/finemapping_merged/ gs://genetics-portal-dev-staging/finemapping/${version_date}_merged


################ PREVIOUS CODE ###############
mkdir -p finemapping_results/190612
# Ed's previous fine-mapping results
gsutil -m rsync -r gs://genetics-portal-staging/finemapping/190612/top_loci finemapping_results/190612/top_loci
mv finemapping_results/190612/top_loci/part-00000-5e9aef09-9aea-4fab-898f-0725f0bbf865-c000.json.gz finemapping_results/190612/top_loci.json.gz

# Newer finemapping results
gsutil -m rsync -r gs://genetics-portal-dev-staging/finemapping finemapping_results/

mkdir finemapping_merged

python merge_finemap_results.py

zcat finemapping_results/190612/top_loci.json.gz \
    finemapping_results/210309/top_loci.json.gz \
    finemapping_results/210515/top_loci.json.gz \
    finemapping_results/finngen_210515/top_loci.json.gz \
    | gzip > finemapping_merged/top_loci.json.gz


echo "Merged fine-mapping results, includes Ed's release + 213 GWAS Catalog + 4 Covid studies + T1D study" > finemapping_merged/README.txt
version_date=`date +%y%m%d`
gsutil -m cp -r finemapping_merged/README.txt gs://genetics-portal-dev-staging/finemapping/merged_$version_date/README.txt
gsutil -m cp -r finemapping_merged/top_loci.json.gz gs://genetics-portal-dev-staging/finemapping/merged_$version_date/top_loci.json.gz
gsutil -m cp -r finemapping_merged/credset gs://genetics-portal-dev-staging/finemapping/merged_$version_date/

```

### Other notes

I did a test run on two different VM instances where I fine-mapped 15 GWAS.
One VM had a balanced persistent disk (200 Gb), one had an SSD (200 Gb). Otherwise they both were N2-standard-8 configurations.
The result was that the SSD version took about 4% longer than the standard disk. I did not try with a local SSD, but I suspect that the disk makes no difference, since the pipeline is CPU-bound.

##### Useful commands

```
# Parse time taken for each run
grep "Time taken" logs/study_id=*/phenotype_id=*/bio_feature=*/chrom=*/logfile.txt
ls -rt logs/study_id=*/phenotype_id=*/bio_feature=*/chrom=*/logfile.txt | xargs grep "Time taken"

# List all
ls logs/study_id=*/phenotype_id=*/bio_feature=*/chrom=*/logfile.txt
```

##### Notes

- Currently fails for sex chromosomes
  - Need to replace X with 23 in plink file or when specifying gcta command
  - Need to impute sex in plink file for X for cojo to work
- Manifest NAs must be represented with "None"
- P-value threshold is specified in 1_scan_input_parquets.py. Set to 5e-8 for GWAS, and (0.05 / num_tests) for mol trait