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

#module load jdk/16.0.1
#module load openjdk/11.0.2
#module load squashfs
#module load singularity

set -f

nextflow_path=[path to Nextflow executable]

eqtls=[Folder with eQTLGen eQTL parquet files]
allele_info=[Allele information in parquet format]

gmt=[Folder with gene set files]
nih_symbols=https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/archive/monthly/tsv/hgnc_complete_set_2023-09-01.txt
# Update this, if neccessary

output_folder=[Output folder]

NXF_VER=23.04.1 ${nextflow_path}/nextflow run HubSnpEnrichmentPipeline.nf \
--eqtls ${eqtls} \
--allele_info ${allele_info} \
--gmt ${gmt} \
--nih_symbols ${nih_symbols} \
--outdir ${output_folder} \
-resume
