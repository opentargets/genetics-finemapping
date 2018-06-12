#

from glob import glob

rule all:
    input:
        "neale_ukb_uk10kref_180502.crediblesets.long.varids.tsv.gz"

rule collate:
    input:
        glob("/nfs/users/nfs_e/em21/otcoregen/em21/finemapping/results/*/credible_set/*.cond.*.credible_sets.tsv")
    output:
        "neale_ukb_uk10kref_180502.crediblesets.long.tsv.gz"
    shell:
        "python scripts/collate_credible_set_locus_files.py --out {output}"

rule annotate_varids:
    """ Maps rsids to variant ids using a neale et al summary stat file
    """
    input:
        credset="neale_ukb_uk10kref_180502.crediblesets.long.tsv.gz",
        sumstat="../uk_biobank_stanford/data/clean/50.nealeUKB_20170915.assoc.clean.tsv.gz"
    output:
        "neale_ukb_uk10kref_180502.crediblesets.long.varids.tsv.gz"
    shell:
        "python scripts/annotate_collated_rsids_w_varids.py "
        "--inf {input.credset} "
        "--outf {output} "
        "--sumstat {input.sumstat}"
