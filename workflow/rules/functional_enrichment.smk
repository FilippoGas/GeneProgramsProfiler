rule run_gsea:
    """
    Run Gene Set Enrichment Analysis on the results of the differential analyses,
    using cytopus gene sets. Produce results for single cell and pseudobulk DEGs
    (both present in the "DEGs_both" csv).
    The outputs for this rule are not defined a priori, as they depend on the number
    of celltypes with DEGs. A .done.txt sentinel is kept as the first output so that
    the R script can locate the output directory; the resulting plots and tables are
    collected through a report directory output.
    """
    input:
        DEGs_both=f"results/{config["analysis_name"]}/DE_analysis/DEGs_both.csv",
        cytopus=f"results/{config["analysis_name"]}/spectra/cytopus_Gene_sets.json",
    output:
        done=f"results/{config["analysis_name"]}/functional_enrichment/run_gsea/.done.txt",
        plots=report(
            directory(
                f"results/{config["analysis_name"]}/functional_enrichment/run_gsea"
            ),
            caption="../report/captions/functional_enrichment/run_gsea.rst",
            category="Functional enrichment",
            subcategory="run_gsea",
            labels={"celltype": "{celltype}"},
            patterns=["{analysis}/{celltype}.pdf", "{analysis}/{celltype}.csv"],
        ),
    log:
        f"logs/{config["analysis_name"]}/functional_enrichment/run_gsea.log",
    benchmark:
        f"results/{config["analysis_name"]}/functional_enrichment/benchmarks/run_gsea.tsv"
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
    The outputs for this rule are not defined a priori, as they depend on the number
    of celltypes with DEGs. A .done.txt sentinel is kept as the first output so that
    the R script can locate the output directory; the resulting plots and tables are
    collected through a report directory output.
    """
    input:
        DEGs_both=f"results/{config["analysis_name"]}/DE_analysis/DEGs_both.csv",
        cytopus=f"results/{config["analysis_name"]}/spectra/cytopus_Gene_sets.json",
        gene_list=f"results/{config["analysis_name"]}/DE_analysis/gene_list.txt",
    output:
        done=f"results/{config["analysis_name"]}/functional_enrichment/run_ora/.done.txt",
        plots=report(
            directory(
                f"results/{config["analysis_name"]}/functional_enrichment/run_ora"
            ),
            caption="../report/captions/functional_enrichment/run_ora.rst",
            category="Functional enrichment",
            subcategory="run_ora",
            labels={"celltype": "{celltype}", "direction": "{direction}"},
            patterns=[
                "{analysis}/{celltype}_{direction}.pdf",
                "{analysis}/{celltype}_{direction}.csv",
            ],
        ),
    log:
        f"logs/{config["analysis_name"]}/functional_enrichment/run_ora.log",
    benchmark:
        f"results/{config["analysis_name"]}/functional_enrichment/benchmarks/run_ora.tsv"
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
