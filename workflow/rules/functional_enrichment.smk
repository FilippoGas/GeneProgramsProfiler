rule run_gsea:
    """
    Run Gene Set Enrichment Analysis on the results of the differential analyses,
    using cytopus gene sets. Produce results for single cell, pseudobulk and 
    merged DEGs.
    """
    input:
    output:
    log:
    conda:
    threads:
    resources:
    params:
    message:
    script: