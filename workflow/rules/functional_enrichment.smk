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
        DEGs_both=f"results/{config["analysis_name"]}/DE_analysis/DEGs_both.csv",
        cytopus=f"results/{config["analysis_name"]}/spectra/cytopus_Gene_sets.json",
    output:
        f"results/{config["analysis_name"]}/functional_enrichment/run_gsea/.done.txt",
    log:
        f"logs/{config["analysis_name"]}/functional_enrichment/run_gsea.log",
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


rule run_ora:
    """
    Run Over Representation Analysis on the results of the differential analyses,
    using cytopus gene sets. Produce results for single cell, pseudobulk and common
    DEGs (all present in the "DEGs_both" csv).
    The output for this rule is not defined a priori, as it depends on the number
    of celltypes with DEGs. For this reason a fake output is used to understand if
    the rule has finished.
    """
    input:
        DEGs_both=f"results/{config["analysis_name"]}/DE_analysis/DEGs_both.csv",
        cytopus=f"results/{config["analysis_name"]}/spectra/cytopus_Gene_sets.json",
        gene_list=f"results/{config["analysis_name"]}/DE_analysis/gene_list.txt",
    output:
        f"results/{config["analysis_name"]}/functional_enrichment/run_ora/.done.txt",
    log:
        f"logs/{config["analysis_name"]}/functional_enrichment/run_ora.log",
    conda:
        "../envs/functional_enrichment.yaml"
    threads: config["functional_enrichment"]["run_ora"]["cores"]
    resources:
        mem_mb=config["functional_enrichment"]["run_ora"]["mem_mb"],
        time=config["functional_enrichment"]["run_ora"]["time"],
        queue=config["queues"]["cpu"],
    params:
        padj_thresh=config["functional_enrichment"]["run_gsea"]["padj_thresh"],
    message:
        "Run Over Representation Analysis on the results of the differential expression analyses."
    script:
        "../scripts/functional_enrichment/run_ora.R"
