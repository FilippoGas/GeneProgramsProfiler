rule make_comp_table:
    """
    Put togheter results from different methods in one table allowing direct comparison

    """
    input:
        cytopus=f"results/{config["analysis_name"]}/spectra/cytopus_Gene_sets.json",
        spectra_GP_activation_WMW=f"results/{config["analysis_name"]}/spectra/spectra_WMW/gp_activation_WMW.csv",
        spectra_GP_activation_LMM=f"results/{config["analysis_name"]}/spectra/spectra_LMM/gp_activation_LMM.csv",
        cNMF_GP_activation_WMW=f"results/{config["analysis_name"]}/cNMF/{config["analysis_name"]}/cNMF_WMW/gp_activation_WMW.csv",
        cNMF_GP_activation_LMM=f"results/{config["analysis_name"]}/cNMF/{config["analysis_name"]}/cNMF_LMM/gp_activation_LMM.csv",
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
        case=config["case_condition"],
        padj_thresh=config["collect_results"]["make_comp_table"]["padj_thresh"],
    message:
        "Put togheter results from different methods in one table allowing direct comparison."
    script:
        "../scripts/collect_results/make_comp_table.R"


rule comp_table_plots:
    """
    Make concordance plots showing detected deregulations across methods.
    """
    input:
        table=f"results/{config["analysis_name"]}/collect_results/comparative_table.csv",
    output:
        horizontal_upset=f"results/{config["analysis_name"]}/collect_results/horizontal_upset.pdf",
        vertical_upset=f"results/{config["analysis_name"]}/collect_results/vertical_upset.pdf",
        intersections=f"results/{config["analysis_name"]}/collect_results/vertical_upset_intersections.xlsx",
    log:
        f"logs/{config["analysis_name"]}/collect_results/comp_table_plots.log",
    conda:
        "../envs/collect_results.yaml"
    threads: config["collect_results"]["comp_table_plots"]["cores"]
    resources:
        mem_mb=config["collect_results"]["comp_table_plots"]["mem_mb"],
        time=config["collect_results"]["comp_table_plots"]["time"],
        queue=config["queues"]["cpu"],
    params:
        FDR_thresh=config["collect_results"]["comp_table_plots"]["FDR_thresh"],
        log2FC_thresh=config["collect_results"]["comp_table_plots"]["log2FC_thresh"],
        effect_size_thresh=config["collect_results"]["comp_table_plots"][
            "effect_size_thresh"
        ],
    message:
        "Make concordance plots showing detected deregulations across methods."
    script:
        "../scripts/collect_results/comp_table_plots.R"
