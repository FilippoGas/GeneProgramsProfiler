rule run_gsea:
    """
    Run Gene Set Enrichment Analysis on the results of the differential analyses,
    using cytopus gene sets. Produce results for single cell and pseudobulk DEGs
    (both present in the "DEGs_both" csv).
    The output for this rule is not defined a priori, as it depends on the number
    of celltypes with DEGs. For this reason a fake output is used to understand if
    the rule has finished.
    """
    input:
        DEGs_both=f"results/DE_analysis/DEGs_both.csv",
        cytopus=f"results/spectra/{config["analysis_name"]}_cytopus_Gene_sets.json",
    output:
        "results/functional_enrichment/run_gsea/.done.txt",
    log:
        "logs/functional_enrichment/run_gsea.log",
    conda:
        "../envs/functional_enrichment.yaml"
    threads: config["functional_enrichment"]["run_gsea"]["cores"]
    resources:
        mem_mb=config["functional_enrichment"]["run_gsea"]["mem_mb"],
        time=config["functional_enrichment"]["run_gsea"]["time"],
        queue=config["queues"]["cpu"],
    params:
        padj_thresh=config["functional_enrichment"]["run_gsea"]["padj_thresh"],
    message:
        "Run Gene Set Enrichment Analysis on the results of the differential expression analyses."
    script:
        "../scripts/functional_enrichment/run_gsea.R"
