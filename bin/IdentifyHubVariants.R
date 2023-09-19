#!/usr/bin/env Rscript
library(argparse)
library(arrow)
library(data.table)
library(dplyr)

setDTthreads(12)

parser <- ArgumentParser(description = 'Extract hub variants.')
parser$add_argument('--sig_eqtls', type = 'character',
                    help = 'File with significant eQTL effects.')
parser$add_argument('--allele_info', type = 'character',
                    help = 'File with SNP information in parquet format.')
parser$add_argument('--gene_count', type = 'numeric', default = 10,
                    help = 'How many genes for hub SNP.')
parser$add_argument('--gene_overlap_thresh', type = 'numeric', default = 0.8,
                    help = 'Prune overlap threshold.')
parser$add_argument('--prune_window', type = 'numeric', default = 1000000,
                    help = 'Prune window.')

args <- parser$parse_args()

# Function prune hub variants
prune_hub_variants <- function(dat, snp_name_col, snp_chr_col, snp_pos_col, gene_col, p_col, hub_gene_th = 10,
                               P_th = 5e-8, window_th = 1000000, overlap_th = 0.8){

  dat <- data.table(SNP = dat[[snp_name_col]],
                        CHR = dat[[snp_chr_col]],
                        POS = dat[[snp_pos_col]],
                        gene = dat[[gene_col]],
                        P = dat[[p_col]])

  message("Data pre-formatted!")

  dat <- dat[P <= P_th]

  dat[, nr_genes := .N, by = SNP]
  dat <- dat[nr_genes >= hub_gene_th]
  setorder(dat, -nr_genes)

  message(paste(length(unique(dat$SNP)), "unpruned hub variants in the input"))

  results <- data.table(SNP = integer(0), CHR = integer(0), POS = integer(0), nr_genes = integer(0))

  while (nrow(dat) > 0) {
    # Get the first SNP
    temp_snp <- dat$SNP[1]
    temp_chr <- dat$CHR[1]
    temp_pos <- dat$POS[1]
    temp_hub <- dat[SNP %in% dat$SNP[1], ]

    # Add to results
    results <- rbind(results, unique(temp_hub[, .(SNP, CHR, POS, nr_genes)]))

    # Remove processed SNP
    dat <- dat[SNP != temp_snp]

    dat_temp <- dat[
      CHR %in% temp_chr & abs(POS - temp_pos) < window_th,
      overlap := sum(gene %in% temp_hub$gene),
      by = .(SNP)
    ]

    dat_temp <- dat_temp[(overlap / nr_genes) < overlap_th]

    dat <- dat[!(CHR %in% temp_chr & abs(POS - temp_pos) < window_th)]
    dat <- unique(rbind(dat, dat_temp))

    # Sort by the number of genes in descending order
    setorder(dat, -nr_genes)
    message(paste(temp_snp, temp_chr, temp_pos, "with", temp_hub$nr_genes[1], "genes", "added!", sep = " "))
  }
  return(results)
}


ZtoP <- function(Z, largeZ = FALSE, log10P = TRUE) {
  if (!is.numeric(Z)) {
    message("Some of the Z-scores are not numbers! Please check why!")
    message("Converting the non-numeric vector to numeric vector.")
    Z <- as.numeric(Z)
  }

  if (largeZ == TRUE) {
    P <- log(2) + pnorm(abs(Z), lower.tail = FALSE, log.p = TRUE)

    if (largeZ == TRUE & log10P == TRUE) {
      P <- -(P * log10(exp(1)))
    }
  } else {
    P <- 2 * pnorm(abs(Z), lower.tail = FALSE)

    if (min(P) == 0) {
      P[P == 0] <- .Machine$double.xmin
      message("Some Z-score indicates very significant effect and P-value is truncated on 2.22e-308. If relevant, consider using largeZ = TRUE argument and logarithmed P-values instead.")
    }
  }

  return(P)
}



# Read in significant effects
message("Read in significant eQTL effects")
eqtls <- fread(args$sig_eqtls, key = "variant")
eqtls$P <- ZtoP(eqtls$beta/eqtls$standard_error)
message("Sig. eQTL effects read in!")
# Read in SNP information
message("Read in SNP information")
alleles <- read_parquet(args$allele_info, col_select = c("ID", "CHR", "bp"))
message("SNP information read in!")
alleles <- setDT(alleles, key = "ID")
message("SNP information converted to data.table!")
# setkey(alleles, "ID")
# message("Keys set!")
alleles <- alleles[ID %in% eqtls$variant]
message("SNP information prefiltered!")
# Merge
eqtls <- merge(eqtls, alleles, by.x = "variant", by.y = "ID")
message("Merged!")
# Find hub variants
message("Find hub variants!")
hub <- prune_hub_variants(eqtls, 
snp_name_col = "variant", 
snp_chr_col = "CHR", 
snp_pos_col = "bp", 
gene_col = "phenotype", 
p_col = "P",
hub_gene_th = args$gene_count,
P_th = max(eqtls$P), 
window_th = 1000000, # TODO add as argument to pipeline options
overlap_th = 0.8) # TODO add as argument to pipeline options
# Write out
if (nrow(hub) == 0){
  message("No hub variants!")
  exit()
}
fwrite(hub, "hub_snps.txt", sep = "\t")
message("Finished!")
