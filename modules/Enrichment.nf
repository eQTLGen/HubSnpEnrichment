#!/bin/bash nextflow

process FILTER {

    //container = 'quay.io/urmovosa/pqtlvseqtl:v0.2'
    
    scratch true
    
    input:
        path(eqtl_folder)
        val(p_thresh)
        val(i2_thresh)

    output:
        path "*_sig_res.txt"

    shell:
        '''
        mkdir tmp_eqtls
        cp -r phenotype=* tmp_eqtls/

        FilterDataBySignificance.R \
        --eqtl_folder tmp_eqtls \
        --i2 !{i2_thresh} \
        --P !{p_thresh}

        rm -r tmp_eqtls
        '''
}

process FINDHUB {

    input:
        tuple path(sig_res), path(alleles), val(gene_overlap), val(gene_overlap_thresh), val(prune_window)

    output:
        path "*.txt"

    script:
        """
        IdentifyHubVariants.R \
        --sig_eqtls ${sig_res} \
        --allele_info ${alleles} \
        --gene_count ${gene_overlap} \
        --gene_overlap_thresh ${gene_overlap_thresh} \
        --prune_window ${prune_window}
        """
}

process EXTRACTHUB {
    
    scratch true

    input:
        path(input)

    output:
        path "*hub_extraction.txt"

    script:
        """
        mkdir tmp_eqtls
        cp -r phenotype=* tmp_eqtls/

        mv *.txt hub_snps2.txt

        ExtractHubVariants.R \
        --eqtl_folder tmp_eqtls \
        --hub_snps hub_snps2.txt

        rm -r tmp_eqtls
        """
}

process SPLIT {

    input:
        path(input)

    output:
        path "*.rds"

    script:
        """
        SplitVariants.R \
        --hub_variants HubVariants.txt
        """
}

process READSETS {

    input:
       tuple path(gmt), path(nih)

    output:
        path "*.rds"

    script:
        """
        PrepareGeneSets.R \
        --gmt ${gmt} \
        --hgnc_map ${nih}
        """
}


process GSEA {

    input:
        tuple path(gmt), path(gene_list)

    output:
        path "*_results.txt"

    script:
        """
        mkdir tmp
        mv *.rds tmp/.
        mv tmp/PrepLibs.rds .

        RunGsea.R \
        --genelist tmp \
        --gmt ${gmt} \
        """
}

workflow FilterDataBySignificance {
    take:
        data

    main:
        FilterDataBySignificance_output_ch = FILTER(data)

    emit:
        FilterDataBySignificance_output_ch
}

workflow IdentifyHubVariants {
    take:
        data

    main:
        IdentifyHubVariants_output_ch = FINDHUB(data)

    emit:
        IdentifyHubVariants_output_ch
}

workflow ExtractHubVariants {
    take:
        data

    main:
        ExtractHubVariants_output_ch = EXTRACTHUB(data)

    emit:
        ExtractHubVariants_output_ch
}

workflow SplitHubs {
    take:
        data

    main:
        Split_output_ch = SPLIT(data)

    emit:
        Split_output_ch
}

workflow PrepareGeneSets {
    take:
        data

    main:
        PrepareGeneSets_output_ch = READSETS(data)

    emit:
        PrepareGeneSets_output_ch
}

workflow RunGsea {
    take:
        data

    main:
        RunGsea_output_ch = GSEA(data)

    emit:
        RunGsea_output_ch
}