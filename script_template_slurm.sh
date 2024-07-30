#!/bin/bash

#SBATCH --time=01:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=6G
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --job-name="HubSnpEnrichment"

# Here load needed system tools (Java 1.8 is required, one of singularity or anaconda - python 2.7 are needed,
# depending on the method for dependancy management)

module load openjdk/11.0.2
module load squashfs
module load singularity

set -f

nextflow_path=[path to nextflow executable]

eqtls=[path to full eQTL mapping files, separated by gene name]
allele_info=[SNP annotation file in parquet format: 1000G-30x.parquet]

gmt=[path to gene libraries used for testing the enrichment]
nih_symbols=https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/archive/monthly/tsv/hgnc_complete_set_2023-09-01.txt

output_folder=[output folder for results]

NXF_VER=23.04.3 ${nextflow_path}/nextflow run HubSnpEnrichmentPipeline.nf \
--eqtls ${eqtls} \
--allele_info ${allele_info} \
--gmt ${gmt} \
--nih_symbols ${nih_symbols} \
--GeneCountThresh 6 \
--GseaBatchSize 1 \
--SnpFilter test_data/test_snp.txt \
--outdir ${output_folder} \
-profile slurm,singularity \
-resume
