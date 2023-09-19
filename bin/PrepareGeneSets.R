#!/usr/bin/env Rscript

library(argparse)
library(dplyr)
library(tidyr)
library(data.table)
library(clusterProfiler)

parser <- ArgumentParser(description = 'Prepare gene sets.')
parser$add_argument('--gmt', type = 'character',
                    help = 'Folder with libraries from enrichr web site.')
parser$add_argument('--hgnc_map', type = 'character',
                    help = 'NIH file for gene annotation.')

args <- parser$parse_args()

# Prepare mapping file
hgnc_map <- fread(args$hgnc_map)
hgnc_map <- hgnc_map[, colnames(hgnc_map) %in% c("ensembl_gene_id", "symbol", "alias_symbol", "prev_symbol"), with = FALSE]
hgnc_map <- hgnc_map %>% unite(symbols, c(symbol, alias_symbol, prev_symbol), sep = "|") %>% 
separate_longer_delim(symbols, delim = "|") %>% filter(ensembl_gene_id != "" & symbols != "")

# read in GMT
message("Read in libraries!")
libs <- list.files(args$gmt, full.names = TRUE)

prep_libs <- list()

for(lib in libs){
message(paste("Preparing:", lib))
gmt <- read.gmt(lib)

libname <- basename(lib)

# Convert HGNC to ENSEMBL
message("Preparing libraries!")

message(paste(table(!unique(gmt$gene) %in% hgnc_map$symbols)["TRUE"], "genes from", length(unique(gmt$gene)), "genes present in gene sets are not in ENSEMBL"))

gmt <- merge(gmt, hgnc_map, by.x = "gene", by.y = "symbols")
gmt <- unique(gmt[, -1])

terms <- unique(gmt$term)
gmt_analysis <- vector(mode = "list", length = length(terms))

for (i in 1:length(terms)){

    temp <- gmt[gmt$term %in% terms[i], ]$ensembl_gene_id
    gmt_analysis[[i]] <- temp

}

names(gmt_analysis) <- paste(libname, unique(gmt$term), sep = "| ")

prep_libs <- c(prep_libs, gmt_analysis)
}

saveRDS(prep_libs, file = "PrepLibs.rds")