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
        table=report(
            f"results/{config["analysis_name"]}/collect_results/comparative_table.csv",
            caption="../report/captions/collect_results/comparative_table.rst",
            category="Methods consensus",
            subcategory="Comparison table",
        ),
    log:
        f"logs/{config["analysis_name"]}/collect_results/make_comp_table.log",
    benchmark:
        f"results/{config["analysis_name"]}/collect_results/benchmarks/make_comp_table.tsv"
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
        horizontal_upset=report(
            f"results/{config["analysis_name"]}/collect_results/horizontal_upset.pdf",
            caption="../report/captions/collect_results/horizontal_upset.rst",
            category="Methods consensus",
            subcategory="comparison plots",
        ),
        vertical_upset=report(
            f"results/{config["analysis_name"]}/collect_results/vertical_upset.pdf",
            caption="../report/captions/collect_results/vertical_upset.rst",
            category="Methods consensus",
            subcategory="comparison plots",
        ),
        intersections=report(
            f"results/{config["analysis_name"]}/collect_results/vertical_upset_intersections.xlsx",
            caption="../report/captions/collect_results/vertical_upset_intersection.rst",
            category="Methods consensus",
            subcategory="comparison plots table",
        ),
    log:
        f"logs/{config["analysis_name"]}/collect_results/comp_table_plots.log",
    benchmark:
        f"results/{config["analysis_name"]}/collect_results/benchmarks/comp_table_plots.tsv"
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


# Benchmark files are not treated as outputs by the report plugin, so they have
# to be exposed explicitly. They are attached as report-flagged inputs of the
# collector rule below, keeping them separated per analysis and per module.
BENCHMARKED_RULES = {
    "preprocess": ["annotate_and_save"],
    "spectra": [
        "prepare_cytopus_list",
        "run_spectra",
        "spectra_rename_programs",
        "spectra_WMW",
        "spectra_WMW_plot",
        "spectra_LMM",
        "spectra_LMM_plot",
    ],
    "DE_analysis": ["run_DE_analysis", "DEA_plots"],
    "functional_enrichment": ["run_gsea", "run_ora"],
    "cNMF": [
        "cNMF_prepare",
        "cNMF_combine",
        "cNMF_k_selection_plot",
        "cNMF_extract_best_k",
        "cNMF_consensus",
        "cNMF_rename_programs",
        "cNMF_WMW",
        "cNMF_WMW_plot",
        "cNMF_LMM",
        "cNMF_LMM_plot",
    ],
    "collect_results": ["make_comp_table", "comp_table_plots"],
}


benchmark_tsv = lambda module, rule_name: report(
    f"results/{config["analysis_name"]}/{module}/benchmarks/{rule_name}.tsv",
    caption="../report/captions/benchmarks/benchmark_table.rst",
    category="Benchmarks",
    subcategory=module,
)


rule collect_benchmarks:
    """
    Collect the benchmarks of all rules into the report, per analysis.
    """
    input:
        benchmarks=[
            *[
                benchmark_tsv(module, rule_name)
                for module, rule_names in BENCHMARKED_RULES.items()
                for rule_name in rule_names
            ],
            *[
                benchmark_tsv("cNMF", f"cNMF_factorize_worker_{worker}")
                for worker in WORKERS
            ],
        ],
    output:
        summary=report(
            f"results/{config["analysis_name"]}/benchmarks/{config["analysis_name"]}_benchmarks.html",
            caption="../report/captions/benchmarks/benchmark_summary.rst",
            category="Benchmarks",
            subcategory="Summary",
        ),
    log:
        f"logs/{config["analysis_name"]}/benchmarks/collect_benchmarks.log",
    conda:
        "../envs/collect_results.yaml"
    threads: config["collect_results"]["collect_benchmarks"]["cores"]
    resources:
        mem_mb=config["collect_results"]["collect_benchmarks"]["mem_mb"],
        time=config["collect_results"]["collect_benchmarks"]["time"],
        queue=config["queues"]["cpu"],
    message:
        "Collect the benchmarks of all rules into the report."
    script:
        "../scripts/collect_results/collect_benchmarks.R"
