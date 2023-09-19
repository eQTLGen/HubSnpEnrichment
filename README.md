# Pipeline for running gene set enrichment analyses (GSEA) for eQTL hub variants

This pipeline runs gene set enrichment analysis (GSEA) for every pleiotropic "hub" eQTL variant, that affects many genes in the genome-wide eQTL analysis. The overall analysis scheme is following:

1. All significant effects are extracted from eQTL data, based on the selected P-value threshold (default: P<5e-8).
2. Potential "hub" variants are extracted based on the number of genes they affect on selected P-value threshold (default: 10).
3. In order to avoid testing many variant proxies, "hub" variants are "pruned" based on affected genes. Starting with the variant with the largest number of affected genes, other variants are removed in the near vicinity (+/-1Mb) if they affect very similar set of genes (by default 80% of genes shared with the "hub" variant). This process is iterated until all very similarly behaving variants are removed from the analysis.
4. All eQTL association summary statistics for each hub variant are extracted from the genome-wide eQTL results.
5. Genes are ordered by eQTL effect strength and direction (eQTL Z-score).
6. GSEA is run for each "hub" SNP.
## Instructions 

TBA

### Requirements

TBA

### Setup of the analysis

TBA

### Input

TBA

### Settings

TBA

### Command

TBA

### Analysis

TBA

### Output

- You need to have Java 8 installed in your HPC environment.
- You need to have Singularity installed in your HPC environment.
- You need SLURM as scheduler.
- You need to download and install Nextflow executable (https://www.nextflow.io/docs/latest/getstarted.html#installation).

## Acknowledgements

This pipeline was written by Urmo Võsa.

Pipeline uses `fgsea` package for running GSEA. `fgsea` is the work of Korotkevich et al.

[Korotkevich, G. et al. Fast gene set enrichment analysis. 060012 Preprint at https://doi.org/10.1101/060012 (2021).](https://www.biorxiv.org/content/10.1101/060012v3)

[Code repo](https://github.com/ctlab/fgsea/)

The original methodology of GSEA:

[Subramanian, A. et al. Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles. Proceedings of the National Academy of Sciences 102, 15545–15550 (2005).](https://www.pnas.org/doi/10.1073/pnas.0506580102)

[Mootha, V. K. et al. PGC-1α-responsive genes involved in oxidative phosphorylation are coordinately downregulated in human diabetes. Nat Genet 34, 267–273 (2003).](https://www.nature.com/articles/ng1180)

