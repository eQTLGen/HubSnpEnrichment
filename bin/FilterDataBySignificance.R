#!/usr/bin/env Rscript

library(arrow)
library(argparse)
library(data.table)
library(dplyr)
library(dbplyr)

# Functions:
PtoZ <- function(Pvals, betas = NULL) {
  Z <- qnorm(Pvals / 2, lower.tail = FALSE)

  if (!is.null(betas)) {

    if (length(Z) != length(betas)){
      stop("Your P-values and betas vectors are of different length!")
    }

    Z[betas < 0] <- -Z[betas < 0]
  }

  return(Z)
}

parser <- ArgumentParser(description = 'Extract hub variants.')
parser$add_argument('--eqtl_folder', type = 'character',
                    help = 'Folder with eQTL results in parquet format.')
parser$add_argument('--i2', type = 'numeric', default = 100,
                    help = 'I2 threshold.')
parser$add_argument('--P', type = 'numeric', default = 5e-8,
                    help = 'Significance P-value threshold.')

args <- parser$parse_args()

Zthresh <- PtoZ(args$P)
I2thresh <- args$i2

ds <- open_dataset(args$eqtl_folder, partitioning = "phenotype", hive_style = TRUE) 
ds <- ds %>% filter(abs(beta/standard_error) >= !!Zthresh & (i_squared <= !!I2thresh | is.na(i_squared)))

ds <- ds %>% collect()

first_gene <- unique(ds$phenotype)[1]

fwrite(ds, paste0(first_gene, "_sig_res.txt"), sep = "\t", quote = FALSE)
