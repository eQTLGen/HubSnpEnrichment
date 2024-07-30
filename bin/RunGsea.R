#!/usr/bin/env Rscript

library(argparse)
library(fgsea)
library(data.table)
library(dplyr)
library(tidyr)

set.seed(2)

parser <- ArgumentParser(description = 'Run GSEA.')
parser$add_argument('--genelist_folder', type = 'character',
                    help = 'Folder with ranked gene lists for hub SNPs.')
parser$add_argument('--gmt', type = 'character',
                    help = 'Gene sets file containing ENSG IDs.')

args <- parser$parse_args()

res_combined <- data.table(
Hub_SNP = NA,
pathway	= NA,
pval = NA,
padj = NA,
log2err = NA,
ES = NA,
NES = NA,
size = NA,
leadingEdge = NA,
GSEA_type = NA,
nr_tested_genes = NA,
nr_pos_dir = NA,
nr_neg_dir = NA,
nr_pos_dir_Z = NA, 
nr_neg_dir_Z = NA)[-1]

fwrite(res_combined, "results.txt", sep = "\t", quote = FALSE)


# Read in prepared gene sets
gmt <- readRDS(args$gmt)

inp_lists <- list.files(args$genelist_folder, full.names = TRUE)

for (i in 1:length(inp_lists)){

if (i == 1) {first_snp <- sub("\\.rds$", "", basename(inp_lists[i]))}

# Hub SNP ID
snp <- sub("\\.rds$", "", basename(inp_lists[i]))

# read in hub SNP genes

genelist <- readRDS(inp_lists[i])

# Run GSEA
message(paste("Running GSEA:", snp))

# Effect size and direction
res <- fgsea(pathways = gmt,
stats = genelist,
minSize  = 20,
maxSize  = 3000)

res <- as.data.table(res)

if (nrow(res) == 0){ res <- data.table(
pathway	= NA,
pval = NA,
padj = NA,
log2err = NA,
ES = NA,
NES = NA,
size = NA,
leadingEdge = NA,
nr_tested_genes = NA,
nr_pos_dir = NA,
nr_neg_dir = NA,
nr_pos_dir_Z = NA, 
nr_neg_dir_Z = NA,
Hub_SNP = NA)}

# Add info about direction
res$GSEA_type <- "Z"
res$nr_tested_genes <- length(genelist)
res$nr_pos_dir <- length(genelist[genelist > 0])
res$nr_neg_dir <- length(genelist[genelist < 0])
res$nr_pos_dir_Z <- length(genelist[genelist > 0 & abs(genelist) > 5.45131])
res$nr_neg_dir_Z <- length(genelist[genelist < 0 & abs(genelist) > 5.45131])
res$Hub_SNP <- snp
# TODO: add possibility to select Z threshold

res <- res[, c(ncol(res), 1:(ncol(res) - 1)), with = FALSE]

res_combined <- rbind(res_combined, res)

# Only effect size
res <- fgsea(pathways = gmt,
stats = abs(genelist),
minSize  = 20,
maxSize  = 3000)

res <- as.data.table(res)

if (nrow(res) == 0){ res <- data.table(
pathway	= NA,
pval = NA,
padj = NA,
log2err = NA,
ES = NA,
NES = NA,
size = NA,
leadingEdge = NA,
nr_tested_genes = NA,
nr_pos_dir = NA,
nr_neg_dir = NA,
nr_pos_dir_Z = NA, 
nr_neg_dir_Z = NA,
Hub_SNP = NA)}

# Add info about direction
res$GSEA_type <- "abs(Z)"
res$nr_tested_genes <- length(genelist)
res$nr_pos_dir <- length(genelist[genelist > 0])
res$nr_neg_dir <- length(genelist[genelist < 0])
res$nr_pos_dir_sig_Z <- length(genelist[genelist > 0 & abs(genelist) > 5.45131])
res$nr_neg_dir_sig_Z <- length(genelist[genelist < 0 & abs(genelist) > 5.45131])
res$Hub_SNP <- snp
# TODO: add possibility to select significant Z threshold

res <- res[, c(ncol(res), 1:(ncol(res) - 1)), with = FALSE]

res_combined <- rbind(res_combined, res)
fwrite(res_combined[, -9, with = FALSE], paste(first_snp, "results.txt", sep = "_"), sep = "\t", quote = FALSE, append = TRUE)
}
