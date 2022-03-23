Convert FinnGen fine-mapping outputs to OT Genetics top loci and credset
==========================================
NOTE: When ingesting a new FinnGen version, you will need to update paths below to refer to the new release.

FinnGen runs fine-mapping using in-sample LD. We could not reproduce a fine-mapping analysis at a similar quality, and so it is better to import their results.
They identify associated regions, and then merge together nearby regions, followed by runnning both FINEMAP and SuSIE. I don't see much reason to include both of these, and we will use the SuSIE output since it is simpler to work with.
FinnGen fine-mapping results can be browsed here:
https://console.cloud.google.com/storage/browser/finngen-public-data-r6/finemapping/

We use only the ".snp.bgz" SuSIE output files, which have the following format:
```
trait	region	v	rsid	chromosome	position	allele1	allele2	maf	beta	se	p	mean	sd	prob	cs
AB1_PEDICULOSIS_ACARIASIS_OTHERINFEST	chr9:96623010-99623010	9:96624087:C:T	chr9_96624087_C_T	chr9	96624087	C	T	0.001635	-0.3375	0.802302	0.674	-9.63943240158685e-08	0.000465118163614926	2.8571
...one line per SNP in the credible set
```

We select the top SNP by p-value within each region for the top_loci table.
We use the finngen_R*_<endpoint>.SUSIE.snp.bgz file for a trait to determine credible sets, since this has all fine-mapped regions for a given trait.

Note: we discussed whether to include a single top variant from each FinnGen fine-mapping region, or to include the top variant from each credible set within a region (as long as the variant had p < 5e-8). We decided to include only the single top variant, because otherwise downstream steps such as colocalisation may not work. In the coloc pipeline (as of 2021-Apr) we use GCTA to get conditionally independent sumstats for each top variant (conditioning on other top_loci variants in the region). With FinnGen we can't do this, since the top variants were determined using SuSIE, and also because we don't have an equivalent reference panel. We wouldn't be able to run colocalisation with conditional sumstats, and so we could only do it using the full signal, which wouldn't correctly represent the secondary/tertiary signals. HOWEVER, in the future FinnGen and eQTL catalogue may both release the direct SuSIE outputs / Bayes Factors for all SNPs (or a filtered subset), and these could be used directly with the new version of coloc to do colocalisations with each independent signal from SuSIE.

### FinnGen R6 stats
2861 endpoints in R6 manifest
2802 studies with fine-mapping outputs (SUSIE .snp.bgz files)
2802 studies with fine-mapping that are in manifest
1291 studies with fine-mapping where credible set is retained (some studies have no GW-sig loci, or the credible set has no SNP that is GW-sig)
6861 top loci

### Requirements
- Spark v2.4.0
- [conda](https://conda.io/docs/)
- GNU parallel

### Setup environment
This step may already have been done for GWAS catalog fine-mapping.
```
git clone https://github.com/opentargets/genetics-finemapping.git
cd ~/genetics-finemapping/finngen

# If fine-mapping pipeline not already set up 
bash setup.sh
conda env create -n finemap --file ../environment.yaml
```

### Configure pipeline

We use the same pipeline config parameters as for the main fine-mapping analysis: `../configs/analysis.config.yaml`

### Import a single study

A single study's SuSIE output can be imported using the function below.

```
# Activate environment
source activate finemap

# View args
$ python finngen_finemapping_ingest.py --help
usage: finngen_finemapping_ingest.py [-h] --file <file>
                               --config_file <file>
                               --study_id <str>
                               --toploci <file>
                               --credset <file>
                               --log <file>

optional arguments:
  -h, --help            show this help message and exit
  --file <file>         Input: finngen_R*_<endpoint>.SUSIE.snp.bgz
  --config_file <file>  Input: analysis config file
  --study_id <str>      study_id corresponding to input file
  --toploci <file>      Output: top loci json file
  --credset <file>      Output: credible set json file
  --log <file>          Output: log file
```

### Running the pipeline
The workflow below assumes a single machine, and runs multiple studies using gnu parallel.

#### Step 1: Prepare input data

```
cd $HOME/genetics-finemapping/finngen
mkdir inputs

# Get list of phenotypes
curl https://r6.finngen.fi/api/phenos | jq -r '.[]| @json' > inputs/r6_finngen.json
# May need to first run: sudo apt install jq

# Copy all SuSIE SNP results into finngen input data folder
mkdir -p $HOME/genetics-finemapping/finngen/data
gsutil -m cp gs://finngen-public-data-r6/finemapping/full/*.SUSIE.snp.bgz $HOME/genetics-finemapping/finngen/data
gsutil -m cp gs://finngen-public-data-r6/finemapping/summaries/*.SUSIE.*.tsv $HOME/genetics-finemapping/finngen/data_summaries

# There will be fewer files than the number of phenotypes, because some
# phenotypes have no significant loci and so no fine-mapping.
ls $HOME/genetics-finemapping/finngen/data/finngen_R6*.SUSIE.snp.bgz > inputs/input_paths_finngen.txt
```

#### Step 2: Make manifest file

The manifest file specifies each study to be imported. The manifest is a JSON lines file with each line containing the following fields:

```json
{
  "study_id": "finngen_R5_AD_LO",
  "in_snp": "~/genetics-finemapping/finngen/data/finemapping/finngen_R5_AD_LO.SUSIE.snp.bgz",
  "out_top_loci": "~/genetics-finemapping/finngen/output/finngen_R5_AD_LO/top_loci.json.gz",
  "out_credset": "~/genetics-finemapping/finngen/output/finngen_R5_AD_LO/credible_set.json.gz",
  "out_log": "~/genetics-finemapping/finngen/output/finngen_R5_AD_LO/logfile.txt",
}
```

The manifest file can be automatically generated using:

```
cd $HOME/genetics-finemapping/finngen

# Creates manifest file, one job per study. Output path `finngen.manifest.json`.
# First edit the script to use the correct input files for the FinnGen version.
python 1_make_finngen_finemapping_manifest.py
```

#### Step 3: Run pipeline

```
mkdir logs

# Test pipeline first
cp finngen.manifest.json finngen.manifest.json.all
cat finngen.manifest.json.all | grep RX_STATIN > finngen.manifest.json
python 2_make_commands.py --quiet
time bash 3_run_commands.sh

tmux   # So run continues if connection is lost

cp finngen.manifest.json.all finngen.manifest.json
# See commands that will be run. Creates commands_todo.txt.gz and commands_done.txt.gz.
python 2_make_commands.py --quiet

# Edit args in `run_commands.sh` (e.g. number of cores) and then
# For R5, took 13 min on 7 cores
NCORES=7
time bash 3_run_commands.sh $NCORES

# Exit tmux with Ctrl+b then d

# Check if all commands have been done
python 2_make_commands.py --quiet
# Some commands may have exited with error, due to there being no SNPs in a credible set.
```

The above command will run all analyses specified in the manifest using GNU parallel. It will create two files `commands_todo.txt.gz` and `commands_done.txt.gz` showing which analyses have not yet/already been done. The pipeline can be stopped at any time and restarted without repeating any completed analyses. You can safely regenerate the `commands_*.txt.gz` commands whilst the pipeline is running using `python 2_make_commands.py --quiet`.

If you get this error:
  ModuleNotFoundError: No module named 'dask'
then I've solved it just by deactivating conda and reactivating. This seems to happen especially when using tmux... I'm not sure why.

#### Step 4: Process the results

```
# Combine the results of all the individual analyses
time python 4_combine_results.py

# If there are commands still to do, check if the same number of loci lacked a credible set
zcat results/logfiles.txt.gz | grep "kept 0 vars in credible" | wc -l
zcat results/logfiles.txt.gz | grep "No SNPs in any credible set" | wc -l
# Sum of lines should be the same as the number of studies in manifest (commands_todo)
zcat commands_todo.txt.gz | wc -l

# Make a note about the run
echo "FinnGen R6 - converted SuSIE credible set output to top_loci and credset formats" > results/README.txt

# Copy the results to GCS
version_date=`date +%y%m%d`
gsutil -m rsync -r $HOME/genetics-finemapping/finngen/results/ gs://genetics-portal-dev-staging/finemapping/finngen_$version_date
```

### (Old) FinnGen R5 stats
2925 endpoints in R5 manifest
2781 studies with fine-mapping outputs (SUSIE .snp.bgz files)
2781 studies with fine-mapping that are in manifest
1267 studies with fine-mapping where credible set is retained (some studies have no GW-sig loci, or the credible set has no SNP that is GW-sig)
5707 top loci

### (Old) FinnGen R4 stats
2264 endpoints in R4 manifest
1115 studies with fine-mapping outputs (SUSIE .snp.bgz files)
981 studies with fine-mapping that are in manifest (strange that some are not in the manifest)
907 studies with fine-mapping where credible set is retained (some GW-sig loci have no credible set listed - I'm not sure why!)

