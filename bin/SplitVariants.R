#!/usr/bin/env Rscript

library(data.table)
library(argparse)

setDTthreads(1)

parser <- ArgumentParser(description = 'Split hub variants.')
parser$add_argument('--hub_variants', type = 'character',
                    help = 'File with hub variants.')

args <- parser$parse_args()

hub <- fread(args$hub_variants)

for (hub_snp in unique(hub$variant)){
temp <- hub[variant %in% hub_snp]
temp <- temp[order(Z, decreasing = TRUE)]

temp_outp <- temp$Z

names(temp_outp) <- temp$phenotype

saveRDS(temp_outp, file = paste0(hub_snp, ".rds"))
}