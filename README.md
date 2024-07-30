# Pipeline for running gene set enrichment analyses (GSEA) for eQTL hub variants

This pipeline runs gene set enrichment analysis (GSEA) for every pleiotropic "hub" eQTL variant, that affects multiple genes in the genome-wide eQTL analysis. For testing the enrichment, the pipeline uses the GSEA implementation from fgsea R package from Korotkevich and colleagues. 

Compared to over-representation analysis where only subset of significantly associated eQTL genes are tested for enrichment, GSEA orders all the tested genes based on eQTL effect size and tests if genes part of the certain gene set (pathway, targets of certain TF, etc.) position higher/lower in the ordered gene list than expected by chance (showing enrichment or depletion). 

The overall analysis scheme is the following:

1. All significant effects are extracted from eQTL data, based on the selected P-value threshold (default: P<5e-8).
2. Potential "hub" variants are extracted based on the number of genes they singificantly affect on selected P-value threshold (default: 10).
3. In order to avoid testing many variant proxies, "hub" variants are "pruned" based on affected genes. Starting with the variant with the largest number of affected genes, other variants are removed in the near vicinity (+/-1Mb), if they affect very similar set of genes (by default 80% of genes shared with the "hub" variant). This process is iterated until all very similarly behaving variants are removed from the analysis.
4. All eQTL association summary statistics for each hub variant are extracted from the genome-wide eQTL results.
5. Genes are ordered by eQTL effect strength and direction (eQTL Z-score). Because it might be informative to also consider direction of eQTL effect, ordering is done two ways: using eQTL Z-score (accounting for eQTL effect direction, i.e. whether hub SNP 'activates' or 'represses' certain gene set or pathway in an uniform way) and using absolute eQTL Z-score (do not consider eQTL effect direction, just whether eQTL genes are part of gene set of pathway).
6. GSEA is run for each "hub" SNP and against all the gene libaries that are input. Multiple testing is currently done per library.
   
## Usage information 

### Requirements of the system

- Have access to HPC with multiple cores.
- Have Bash >=3.2 installed.
- Have Java >=11 installed.
- Have Slurm scheduler managing the jobs in the HPC.
- HPC has Singularity installed and running.

### Setup of the pipeline

You can either clone it by using git (if available in HPC):

git clone https://github.com/eQTLGen/HubSnpEnrichment.git

Or just download this from the github download link and unzip.

#### Input

`--eqtls`                       Folder of genome-wide eQTL meta-analysis parquet files. Output of `MetaAnalysis` pipeline.

`--allele_info`                 Parquet file with variant allele info for eQTLGen dataset.

`--gmt`                         Folder with gene GMT (gene matrix transposed) format gene libraries for which to test the enrichment for. Can contain one or multiple GMT files.

`--nih_symbols`                 File with gene symbol aliases and corresponding ENSEMBL gene IDs.

#### Settings

`--outdir`                    Output directory. Defaults to "results".

`--Pthresh`                   P-value threshold for eQTL SNP significance. Defaults to 5e-8.

`--I2thresh `                 Heterogeneity I2 threshold for eQTL inclusion. Defaults to 40.

`--GeneCountThresh `          Number of eGenes to consider variant be "hub". Defaults to 10.

`--GeneOverlapThresh  `       Overlap threshold for pruning hub variants based on eGene overlap. Defaults to 0.8.

`--PruneWindow`               Genomic window used for pruning hub variants based on eGene overlap. Defaults to 1000000.

`--GseaBatchSize`             Argument specifing how many hub SNPs to run through GSEA in a loop. Larger number parallelizes less but increases  Defaults to 5.

`--SnpFilter`                 SNP filtering file, useful for confining the analysis for preselected set of variants. 

#### Command

```sh
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
```

### Output

The pipeline writes out the file `GSEA_results.txt` that contains following columns:

- `Hub_SNP`           Name of the hub variant.
- `pathway`           The pathway/gene set that was tested.
- `pval`              P-value from GSEA.
- `padj`              P-value adjusted for multiple testing by Benjamini-Hochberg FDR. NB! Currently it is adjusted per each library file. If you include multiple GMT files into the analysis you might need to recalculate it.
- `log2err`           The expected error for the standard deviation of the P-value logarithm, as per fgwas documentation.
- `ES`                Enrichment score from GSEA.
- `NES`               Enrichment score normalised to mean enrichment of random samples of the same size, as per fgwas documentation.
- `size`              Size of the pathway after removing gene not present in eQTL dataset.
- `GSEA_type`         How was eGene ordering done (Z or abs(Z)).
- `nr_tested_genes`   How many genes were present in eQTL dataset. 
- `nr_pos_dir`        How many genes showed increased eQTL effect for eQTL effect allele. 
- `nr_neg_dir`        How many genes showed decreased eQTL effect for eQTL effect allele. 
- `nr_pos_dir_sig_Z`  How many genes showed increased eQTL effect for eQTL effect allele AND absolute Z>5.451. 
- `nr_neg_dir_sig_Z`  How many genes showed decreased eQTL effect for eQTL effect allele AND absolute Z>5.451.

**Note:** One of the enrichment analyses could be to test enrichment for transcription factor targets. For such analysis it is informative to test if corresponding transcription factor is also encoded near the corresponding hub SNP, pointing to potential mechanism causing eQTL effect(s). To do so, there is helper script that can be run by `run_annotate_script.sh`. This script writes out `GSEA_results_TF_enrichment_annotated.txt.gz` that contains additional columns with TF information and column `TF_close` indicating whether TF whose targets are enriched among hub SNP eQTL genes is encoded near the corresponding hub SNP (no/yes). The window between hub variant and TF defaults to 1Mb but can be specified in the settings. Script assumes that gene set names in the library start with TF name in HGNC gene symbol format (e.g "TP53 ENCODE").

## Acknowledgements

This pipeline was written by Urmo Võsa.

Pipeline uses `fgsea` package for running GSEA. `fgsea` is the work of Korotkevich et al.

[Korotkevich, G. et al. Fast gene set enrichment analysis. 060012 Preprint at https://doi.org/10.1101/060012 (2021).](https://www.biorxiv.org/content/10.1101/060012v3)

[Code repo](https://github.com/ctlab/fgsea/)

The original methodology of GSEA:

[Subramanian, A. et al. Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles. Proceedings of the National Academy of Sciences 102, 15545–15550 (2005).](https://www.pnas.org/doi/10.1073/pnas.0506580102)

[Mootha, V. K. et al. PGC-1α-responsive genes involved in oxidative phosphorylation are coordinately downregulated in human diabetes. Nat Genet 34, 267–273 (2003).](https://www.nature.com/articles/ng1180)

