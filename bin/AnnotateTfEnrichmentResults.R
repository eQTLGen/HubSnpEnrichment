library(data.table)
library(tidyverse)
library(rtracklayer)
library(arrow)
library(tidyr)
library(argparse)

parser <- ArgumentParser(description = 'Annotate Hub SNP TF enrichment results. Adds information whether SNP is near the corresponding TF.')
parser$add_argument('--gsea_result_file', type = 'character',
                    help = 'File with GSEA results. Can be gzipped.')
parser$add_argument('--hgnc_map', type = 'character',
                    help = 'File with ENSG, HGNC gene names and their aliases. E.g. from here: https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/archive/monthly/tsv/hgnc_complete_set_2023-09-01.txt')
parser$add_argument('--gmt', type = 'character',
                    help = 'Ensembl GMT file with gene coordinates.')
parser$add_argument('--snp_ref', type = 'character',
                    help = 'Variant reference file used in eQTLGen p2 project.')
parser$add_argument('--snp_tf_window', type = "integer", default = 1000000, 
                    help = 'Distance, how far the TSS of TF needs to be from hub SNP to consider it to be closeby.')

args <- parser$parse_args()

tf <- fread(args$gsea_result_file)
#tf <- tf[, -9, with = FALSE]

# Gene aliases
hgnc_map <- fread(args$hgnc_map)
hgnc_map <- hgnc_map[, colnames(hgnc_map) %in% c("ensembl_gene_id", "symbol", "alias_symbol", "prev_symbol"), with = FALSE]
hgnc_map <- hgnc_map %>% unite(symbols, c(symbol, alias_symbol, prev_symbol), sep = "|") %>%
  separate_longer_delim(symbols, delim = "|") %>% filter(ensembl_gene_id != "" & symbols != "")

# Replace TF with ENSG:

TF <- str_replace(tf$pathway, ".*\\| ", "")
TF <- str_replace(TF, " .*", "")

tf$TF <- TF

gmt <- readGFF(args$gmt)
gmt <- gmt[gmt$type == "gene", ]
gmt <- gmt[, c(1, 4, 5, 9)]
gmt <- unique(gmt)

tf <- merge(tf, hgnc_map, by.x = "TF", by.y = "symbols")
tf <- merge(tf, gmt, by.x = "ensembl_gene_id", by.y = "gene_id")

colnames(tf)[17:19] <- c("TF_chr", "TF_start", "TF_end")
# SNP coordinates

snp <- open_dataset(args$snp_ref)
snp_coord <- snp %>% filter(ID %in% tf$Hub_SNP) %>% select(c(ID, bp, CHR)) %>% collect()

tf <- merge(tf, snp_coord, by.x = "Hub_SNP", by.y = "ID")
colnames(tf)[c(20, 21)] <- c("SNP_POS", "SNP_CHR")

tf$TF_close <- "no"
tf[SNP_CHR == TF_chr & (abs(TF_start - SNP_POS) < args$snp_tf_window | abs(TF_end - SNP_POS) < args$snp_tf_window)]$TF_close <- "yes"

outp_file_name <- str_replace(args$gsea_result_file, "\\.gz", "") %>% str_replace("\\.txt", "")
outp_file_name <- paste0(outp_file_name, "_TF_enrichment_annotated.txt.gz")

message("Writing results...")
message(paste("File name:", outp_file_name))
fwrite(tf, file = outp_file_name, sep = "\t")
message("Writing results...done!")
