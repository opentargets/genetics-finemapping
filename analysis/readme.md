Open Target Genetics fine-mapping pipeline
==========================================

The code here is for running fine-mapping analyses using FINEMAP for a handfule of traits to compare against the results obtained using the normal pipeline (GCTA + ABF fine-mapping).
It assumes that you have set up the environment as for our normal fine-mapping pipeline (see readme.md in parent folder).

```
# Activate environment
source activate finemap
```

### Step 1: Copy down LD and GWAS data for fine-mapping
```
cd genetics-finemapping

# Initially use 1000 genomes due to UKB security issues
mkdir -p /home/js29/genetics-finemapping/data/1000Genomes_phase3/EUR/
gsutil -m cp -r gs://genetics-portal-input/1000Genomes_phase3/plink_format_b38/EUR/*  /home/js29/genetics-finemapping/data/1000Genomes_phase3/EUR/

# UKB LD
gsutil -m cp /lustre/scratch115/realdata/mdt3/projects/otcoregen2/em21/uk_biobank_analysis/create_10k_subsample/output/ukb_v3_downsampled10k_plink/* gs://genetics-portal-analysis/js29/ukb/tmp/
gsutil rm gs://genetics-portal-analysis/js29/ukb/tmp/*

mkdir -p /home/js29/genetics-finemapping/data/ukb_downsampled10k/
gsutil -m cp -r gs://genetics-portal-analysis/js29/ukb/tmp/*  /home/js29/genetics-finemapping/data/ukb_downsampled10k/

# A handful of traits were run through the significant window extraction, which
# was changed to also save the intervals as a separate TSV file.
mkdir -p /home/js29/genetics-finemapping/data/filtered/significant_window_2mb_analysis/gwas/
gsutil -m cp -r gs://genetics-portal-sumstats-b38/filtered/significant_window_2mb_analysis/gwas/*  /home/js29/genetics-finemapping/data/filtered/significant_window_2mb_analysis/gwas/
```

### Step 2: Make manifest file

First manually edit 2_make_manifest.py to use the 1000 genomes LD files.

```
# Scan the input summary statistics
python 1_scan_input_parquets.py

# Create manifest file. Output path `configs/manifest.json.gz`
python 2_make_manifest.py
```
### Step 3: Run pipeline

```
# Edit args in `4_run_commands.sh` then
tmux
bash 4_run_commands.sh
```

### Step 4: Remake manifest file for UKB

Manually edit 2_make_manifest.py to use the UKB LD files.

```
# Create manifest file. Output path `configs/manifest.json.gz`
python 2_make_manifest.py
```

### Step 5: Run pipeline again

```
bash 4_run_commands.sh
```

### Step 6: Compare GCTA and FINEMAP credible sets

```
Rscript finemap_gcta_comparison.Rmd
```

