import os


rule make_comp_table:
    """
    Put togheter results from different methods in one table allowing direct comparison

    """
    input:
        cytopus=f"results/{config["analysis_name"]}/spectra/cytopus_Gene_sets.json",
        spectra_GP_activation_WMW=f"results/{config["analysis_name"]}/spectra/spectra_WMW/gp_activation_WMW.csv",
        spectra_GP_activation_LMM=f"results/{config["analysis_name"]}/spectra/spectra_LMM/gp_activation_LMM.csv",
        cNMF_GP_activation_WMW=f"results/{config["analysis_name"]}/cNMF/cNMF_WMW/gp_activation_WMW.csv",
        cNMF_GP_activation_LMM=f"results/{config["analysis_name"]}/cNMF/cNMF_LMM/gp_activation_LMM.csv",
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
        enrichment_dir=lambda wildcards, input: input.gsea.split("run_gsea")[0],
    message:
        "Put togheter results from different methods in one table allowing direct comparison."
    script:
        "../scripts/collect_results/make_comp_table.R"
