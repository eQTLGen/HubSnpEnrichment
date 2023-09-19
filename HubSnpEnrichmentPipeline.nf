#!/usr/bin/env nextflow

/*
 * enables modules
 */
nextflow.enable.dsl = 2
def helpmessage() {

log.info"""

HubSnpEnrichment v${workflow.manifest.version}"
===========================================================
Pipeline for running gene set enrichment analyses for hub SNPs.

Usage:

nextflow run HubSnpEnrichmentPipeline.nf \
--eqtls \
--allele_info \
--gmt \
--outdir

Mandatory arguments:
--eqtls                     eQTLGen parquet dataset.
--allele_info               File with allele info for eQTLGen dataset.
--gmt                       GMT link with gene sets to test.
--nih_symbols               ENSEMBL GTF file for gene annotation.

Optional arguments:
--outdir                    Output directory. Defaults to "results".
--Pthresh                   P-value threshold for eQTL SNP significance. Defaults to 5e-8.
--I2thresh                  Heterogeneity I2 threshold for eQTL inclusion. Defaults to 40.
--GeneCountThresh           Number of eGenes to consider variant be "hub". Defaults to 10.
--GeneOverlapThresh         Overlap threshold for pruning hub variants based on eGene overlap. Defaults to 0.8.
--PruneWindow               Genomic window used for pruning hub variants based on eGene overlap. Defaults to 1000000.
--GseaBatchSize             Argument specifing how many hub SNPs to run through GSEA in a loop. Larger number parallelizes less but increases  Defaults to 5.

""".stripIndent()

}

if (params.help){
    helpmessage()
    exit 0
}

// Default parameters
params.outdir = 'results'
params.inclusion_step_output = 'NO_FILE'
params.Pthresh = 5e-8
params.I2thresh = 40
params.GeneCountThresh = 10
params.GeneOverlapThresh = 0.8
params.PruneWindow = 1000000
params.GseaBatchSize = 5

allele_ch = Channel.fromPath(params.allele_info)
gmt_ch = Channel.fromPath(params.gmt)
nih_ch = Channel.fromPath(params.nih_symbols)

// Define parameter channels
p_thresh = Channel.value(params.Pthresh)
i2_thresh = Channel.value(params.I2thresh)
gene_count_thresh = Channel.value(params.GeneCountThresh)
gene_overlap_thresh = Channel.value(params.GeneOverlapThresh)
prune_window = Channel.value(params.PruneWindow)

//Show parameter values
log.info """=======================================================
HubSnpEnrichment v${workflow.manifest.version}"
======================================================="""
def summary = [:]
summary['Pipeline Version']                         = workflow.manifest.version
summary['Current user']                             = "$USER"
summary['Current home']                             = "$HOME"
summary['Current path']                             = "$PWD"
summary['Working dir']                              = workflow.workDir
summary['Script dir']                               = workflow.projectDir
summary['Config Profile']                           = workflow.profile
summary['Container Engine']                         = workflow.containerEngine
if(workflow.containerEngine) summary['Container']   = workflow.container
summary['Output directory']                         = params.outdir
summary['eQTL results']                             = params.eqtls
summary['Gene set file']                            = params.gmt
summary['Symbol linking']                           = params.nih_symbols
summary['Allele info file']                         = params.allele_info
summary['P-value threshold']                        = params.Pthresh
summary['I2 threshold']                             = params.I2thresh
summary['Gene count threshold']                     = params.GeneCountThresh
summary['Gene overlap threshold']                   = params.GeneOverlapThresh
summary['Gene prune window']                        = params.PruneWindow
summary['Batch size']                               = params.GseaBatchSize

// import modules
include { FILTER; FINDHUB; EXTRACTHUB; GSEA; SPLIT; READSETS; FilterDataBySignificance; IdentifyHubVariants; ExtractHubVariants; SplitHubs; RunGsea; PrepareGeneSets } from './modules/Enrichment.nf'

log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "======================================================="

// Define argument channels
// Get eQTL channel
eqtls = Channel.fromPath(params.eqtls + '/phenotype=*', type: 'dir').buffer(size: 24, remainder: true)

workflow {

        FILTER(eqtls, p_thresh, i2_thresh)
        findhub_ch = FILTER.out.flatten().collectFile(name: 'SigEffects.txt', keepHeader: true, sort: true, storeDir: "${params.outdir}")
        .combine(allele_ch)
        .combine(gene_count_thresh)
        .combine(gene_overlap_thresh)
        .combine(prune_window)
        FINDHUB(findhub_ch)
        EXTRACTHUB(eqtls.combine(FINDHUB.out))
        SPLIT(EXTRACTHUB.out.flatten().collectFile(name: 'HubVariants.txt', keepHeader: true, sort: true))
        READSETS(gmt_ch.combine(nih_ch))
        GSEA(READSETS.out.combine(SPLIT.out.flatten().buffer(size: params.GseaBatchSize, remainder: true)))
        GSEA.out.flatten().collectFile(name: 'GSEA_results.txt', keepHeader: true, sort: true, storeDir: "${params.outdir}")
}
workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
