#!/usr/bin/env Rscript

library(arrow)
library(argparse)
library(data.table)
library(dplyr)

parser <- ArgumentParser(description = 'Extract hub variants.')
parser$add_argument('--eqtl_folder', type = 'character',
                    help = 'Folder with eQTL results.')
parser$add_argument('--hub_snps', type = 'character',
                    help = 'File with hub SNP rs IDs.')

args <- parser$parse_args()

hub_snps <- fread(args$hub_snps)
hub_snps <- hub_snps$SNP

ds <- arrow::open_dataset(args$eqtl_folder, partitioning = "phenotype", hive_style = TRUE)

hubs <- ds %>% filter(variant_index %in% !!hub_snps) %>% collect() %>% as.data.table()

hubs$Z <- hubs$beta / hubs$standard_error
hubs <- hubs[, colnames(hubs) %in% c("variant_index", "phenotype", "Z"), with = FALSE]

hubs <- hubs[order(hubs$Z, decreasing = TRUE)]

fwrite(hubs, paste0(unique(hubs$variant_index)[1], "hub_extraction.txt"), sep = "\t")
