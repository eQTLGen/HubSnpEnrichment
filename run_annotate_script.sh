Rscript --vanilla bin/AnnotateTfEnrichmentResults.R \
--gsea_result_file [Result file of the main pipeline] \
--hgnc_map [Mapping file with HGNC gene name synonyms] \
--gmt [ENSEMBL .gmt file] \
--snp_ref [Variant reference file in .parquet format] \
--snp_tf_window 1000000 