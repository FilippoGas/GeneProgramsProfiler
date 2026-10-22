rule prepare_cytopus_list:
    """
    Download genesets from cytopus for the cell types present in the scRNAseq dataset.
    For more details visit cytopus (https://github.com/wallet-maker/cytopus)
    """
    input:
        config["celltype_conversion_dictionary"],
    output:
        f"results/{config["analysis_name"]}/spectra/cytopus_Gene_sets.json",  # FIXME: refactor output paths so analysis_name sits under results/
    log:
        f"logs/{config["analysis_name"]}/spectra/prepare_cytopus_list.log",
    benchmark:
        f"results/{config["analysis_name"]}/spectra/benchmarks/prepare_cytopus_list.tsv"
    conda:
        "../envs/spectra.yaml"
    threads: config["spectra"]["prepare_cytopus_list"]["cores"]
    resources:
        mem_mb=config["spectra"]["prepare_cytopus_list"]["mem_mb"],
        time=config["spectra"]["prepare_cytopus_list"]["time"],
        queue=config["queues"]["cpu"],
    params:
        global_celltype=config["spectra"]["prepare_cytopus_list"]["global_celltype"],
    message:
        "Retrieve cytopus gene sets for the cell types in the dataset."
    script:
        "../scripts/spectra/prepare_cytopus_list.py"


rule run_spectra:
    """
    Run spectra to quantify gene programs activation in the scRNAseq dataset
    """
    input:
        sc_dataset=f"results/{config["analysis_name"]}/preprocess/AnnData/adata.h5ad",
        cytopus_list=f"results/{config["analysis_name"]}/spectra/cytopus_Gene_sets.json",
    output:
        gene_scores=f"results/{config["analysis_name"]}/spectra/spectra_output/gene_scores_labmda_{config["spectra"]["run_spectra"]["lambda"]}.csv",
        cell_scores=f"results/{config["analysis_name"]}/spectra/spectra_output/cell_scores_labmda_{config["spectra"]["run_spectra"]["lambda"]}.csv",
        factor_markers=f"results/{config["analysis_name"]}/spectra/spectra_output/factor_markers_labmda_{config["spectra"]["run_spectra"]["lambda"]}.csv",
    log:
        f"logs/{config["analysis_name"]}/spectra/run_spectra.log",
    benchmark:
        f"results/{config["analysis_name"]}/spectra/benchmarks/run_spectra.tsv"
    conda:
        "../envs/spectra.yaml"
    threads: config["spectra"]["run_spectra"]["cores"]
    resources:
        mem_mb=config["spectra"]["run_spectra"]["mem_mb"],
        time=config["spectra"]["run_spectra"]["time"],
        queue=config["queues"]["cpu"],
    params:
        lam=config["spectra"]["run_spectra"]["lambda"],
    message:
        "Run spectra to identify activated gene programs in the single cell RNA seq dataset."
    script:
        "../scripts/spectra/run_spectra.py"


rule spectra_rename_programs:
    """
    Spectra might fail to label some factors (gene programs) due to low overlap with
    cytopus gene sets. Here we try to label them through ORA enrichment of the marker
    genes of each unlabeled factor
    """
    input:
        cell_scores=f"results/{config["analysis_name"]}/spectra/spectra_output/cell_scores_labmda_{config["spectra"]["run_spectra"]["lambda"]}.csv",
        factor_markers=f"results/{config["analysis_name"]}/spectra/spectra_output/factor_markers_labmda_{config["spectra"]["run_spectra"]["lambda"]}.csv",
        cytopus_list=f"results/{config["analysis_name"]}/spectra/cytopus_Gene_sets.json",
        metadata=f"results/{config["analysis_name"]}/preprocess/metadata.csv",
    output:
        cell_scores=f"results/{config["analysis_name"]}/spectra/spectra_output/cell_scores_labmda_{config["spectra"]["run_spectra"]["lambda"]}_known_programs.csv",
    log:
        f"logs/{config["analysis_name"]}/spectra/rename_programs.log",
    benchmark:
        f"results/{config["analysis_name"]}/spectra/benchmarks/spectra_rename_programs.tsv"
    conda:
        "../envs/spectra.yaml"
    threads: config["spectra"]["rename_programs"]["cores"]
    resources:
        mem_mb=config["spectra"]["rename_programs"]["mem_mb"],
        time=config["spectra"]["rename_programs"]["time"],
        queue=config["queues"]["cpu"],
    message:
        "Label unknown programs in the cell scores matrix."
    script:
        "../scripts/spectra/rename_programs.R"


rule spectra_WMW:
    """
    Analyze spectra results and test for differential activation of gene programs
    between the two conditions. Gene programs activation differences are
    tested with Wilcoxon-Mann-Whitney U-test.
    """
    input:
        cell_scores=f"results/{config["analysis_name"]}/spectra/spectra_output/cell_scores_labmda_{config["spectra"]["run_spectra"]["lambda"]}_known_programs.csv",
    output:
        gp_activation_WMW=f"results/{config["analysis_name"]}/spectra/spectra_WMW/gp_activation_WMW.csv",
    log:
        f"logs/{config["analysis_name"]}/spectra/spectra_WMW.log",
    benchmark:
        f"results/{config["analysis_name"]}/spectra/benchmarks/spectra_WMW.tsv"
    conda:
        "../envs/spectra.yaml"
    threads: config["spectra"]["spectra_WMW"]["cores"]
    resources:
        mem_mb=config["spectra"]["spectra_WMW"]["mem_mb"],
        time=config["spectra"]["spectra_WMW"]["time"],
        queue=config["queues"]["cpu"],
    params:
        control=config["control_condition"],
        case=config["case_condition"],
        active_cell_thresh=config["spectra"]["spectra_WMW"]["active_cell_thresh"],
    message:
        "Analyze spectra's results and test for differential program activation with Wilcoxon-Mann-Whitney U-test."
    script:
        "../scripts/differential_activation/WMW.R"


rule spectra_WMW_plot:
    """
    Plot results of the differential activation analysis of spectra's output.
    The number of outputs for this rule is not defined a priori, as it depends on the
    number of deregulated programs and of celltypes with DEGs. A .done.txt sentinel is
    kept as the first output so that the R script can locate the output directory; the
    resulting plots are collected through a report directory output.
    """
    input:
        cell_scores=f"results/{config["analysis_name"]}/spectra/spectra_output/cell_scores_labmda_{config["spectra"]["run_spectra"]["lambda"]}_known_programs.csv",
        gp_activation_WMW=f"results/{config["analysis_name"]}/spectra/spectra_WMW/gp_activation_WMW.csv",
    output:
        done=f"results/{config["analysis_name"]}/spectra/spectra_WMW_plots/.done.txt",
        plots=report(
            directory(f"results/{config["analysis_name"]}/spectra/spectra_WMW_plots"),
            caption="../report/captions/spectra/spectra_WMW_plot.rst",
            category="Spectra",
            subcategory="spectra_WMW_plot",
            labels={"celltype": "{celltype}"},
            patterns=["{celltype}.pdf"],
        ),
    log:
        f"logs/{config["analysis_name"]}/spectra/spectra_WMW_plots.log",
    benchmark:
        f"results/{config["analysis_name"]}/spectra/benchmarks/spectra_WMW_plot.tsv"
    conda:
        "../envs/spectra.yaml"
    threads: config["spectra"]["spectra_WMW_plots"]["cores"]
    resources:
        mem_mb=config["spectra"]["spectra_WMW_plots"]["mem_mb"],
        time=config["spectra"]["spectra_WMW_plots"]["time"],
        queue=config["queues"]["cpu"],
    params:
        effect_size_thresh=config["spectra"]["spectra_WMW_plots"]["effect_size_thresh"],
        FDR_thresh=config["spectra"]["spectra_WMW_plots"]["FDR_thresh"],
    message:
        "Plot results from differential gene program activation analysis (WMW)."
    script:
        "../scripts/differential_activation/WMW_plots.R"


rule spectra_LMM:
    """
    Analyze spectra results and test for differential activation of gene programs
    between the two conditions. Gene programs activation differences are
    tested with Linear Mixed Models to correct for cell cycle phase.
    """
    input:
        cell_scores=f"results/{config["analysis_name"]}/spectra/spectra_output/cell_scores_labmda_{config["spectra"]["run_spectra"]["lambda"]}_known_programs.csv",
    output:
        gp_activation_LMM=f"results/{config["analysis_name"]}/spectra/spectra_LMM/gp_activation_LMM.csv",
    log:
        f"logs/{config["analysis_name"]}/spectra/spectra_LMM.log",
    benchmark:
        f"results/{config["analysis_name"]}/spectra/benchmarks/spectra_LMM.tsv"
    conda:
        "../envs/spectra.yaml"
    threads: config["spectra"]["spectra_LMM"]["cores"]
    resources:
        mem_mb=config["spectra"]["spectra_LMM"]["mem_mb"],
        time=config["spectra"]["spectra_LMM"]["time"],
        queue=config["queues"]["cpu"],
    params:
        control=config["control_condition"],
        case=config["case_condition"],
        active_cell_thresh=config["spectra"]["spectra_LMM"]["active_cell_thresh"],
    message:
        "Analyze spectra's results and test for differential program activation with Linear Mixed Models."
    script:
        "../scripts/differential_activation/LMM.R"


rule spectra_LMM_plot:
    """
    Plot results of the differential activation analysis of spectra's output, testing
    differential activation with linear mixed models corrected for cell cycle phase.
    The volcano plot is the only output of this rule and is declared explicitly.
    """
    input:
        gp_activation_LMM=f"results/{config["analysis_name"]}/spectra/spectra_LMM/gp_activation_LMM.csv",
    output:
        volcano=report(
            f"results/{config["analysis_name"]}/spectra/spectra_LMM_plots/GP_deregulation_LMM_volcano.pdf",
            caption="../report/captions/spectra/spectra_LMM_plot.rst",
            category="Spectra",
            subcategory="spectra_LMM_plot",
        ),
    log:
        f"logs/{config["analysis_name"]}/spectra/spectra_LMM_plots.log",
    benchmark:
        f"results/{config["analysis_name"]}/spectra/benchmarks/spectra_LMM_plot.tsv"
    conda:
        "../envs/spectra.yaml"
    threads: config["spectra"]["spectra_LMM_plots"]["cores"]
    resources:
        mem_mb=config["spectra"]["spectra_LMM_plots"]["mem_mb"],
        time=config["spectra"]["spectra_LMM_plots"]["time"],
        queue=config["queues"]["cpu"],
    params:
        log2FC_thresh=config["spectra"]["spectra_LMM_plots"]["log2FC_thresh"],
        FDR_thresh=config["spectra"]["spectra_LMM_plots"]["FDR_thresh"],
    message:
        "Plot results from differential gene program activation analyis. (LMM)"
    script:
        "../scripts/differential_activation/LMM_plots.R"
