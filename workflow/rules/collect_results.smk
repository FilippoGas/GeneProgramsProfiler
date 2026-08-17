import os

rule make_comp_table:
    """
    Put togheter results from different methods in one table allowing
    direct comparison between different methods
    """
    input:
        cytopus=f"results/{config["analysis_name"]}/spectra/cytopus_Gene_sets.json",
        gp_activation_WMW=f"results/{config["analysis_name"]}/spectra/spectra_WMW/gp_activation_WMW.csv",
        gsea=f"results/{config["analysis_name"]}/functional_enrichment/run_gsea/.done.txt",
        ora=f"results/{config["analysis_name"]}/functional_enrichment/run_ora/.done.txt",
    output:
        table=f"results/{config["analysis_name"]}/collect_results/comparative_table.csv",
    log:
        f"logs/{config["analysis_name"]}/collect_results/make_comp_table.log",
    conda:
        "../envs/collect_results.yaml"
    threads: config["collect_results"]["make_comp_table"]["cores"]
    resources:
        mem_mb=config["collect_results"]["make_comp_table"]["mem_mb"],
        time=config["collect_results"]["make_comp_table"]["time"],
        queue=config["queues"]["cpu"],
    params:
        enrichment_dir=f"results/{config["analysis_name"]}/functional_enrichment/
    script:
        "../scripts/collect_results/make_comp_table.R"
