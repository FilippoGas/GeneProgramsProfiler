rule run_DE_analysis:
    """
    Perform single cell and pseudobulk differential expression analysis on the
    input scRNAseq dataset and save DEGs as csv. Save a csv for each analysis type
    (sc and pseudobulk) and one for genes that are deregulated according to both.
    """
    input:
        sc_dataset=f"results/preprocess/rds/{config["analysis_name"]}.rds",
    output:
        DEGs_sc=f"results/DE_analysis/DEGs_sc.csv",
        DEGs_pb=f"results/DE_analysis/DEGs_pb.csv",
        DEGs_both=f"results/DE_analysis/DEGs_both.csv",
    log:
        "logs/DE_analysis/run_DE_analysis.log",
    conda:
        "../envs/DE_analysis.yaml"
    threads: config["DE_analysis"]["run_DE_analysis"]["cores"]
    resources:
        mem_mb=config["DE_analysis"]["run_DE_analysis"]["mem_mb"],
        time=config["DE_analysis"]["run_DE_analysis"]["time"],
        queue=config["queues"]["cpu"],
    params:
        logFC_thresh=config["DE_analysis"]["run_DE_analysis"]["logFC"],
        FDR_thresh=config["DE_analysis"]["run_DE_analysis"]["FDR"],
        control=config["control_condition"],
        case=config["case_condition"],
        sample_col=config["DE_analysis"]["run_DE_analysis"]["sample_column"],
        condition_col=config["DE_analysis"]["run_DE_analysis"]["condition_column"],
    message:
        "Run differential expression analysis in single cell and pseudo bulk mode on the scRNAseq dataset."
    script:
        "../scripts/DE_analysis/run_DE_analysis.R"


# TODO improve plots quality
rule DEA_plots:
    """
    Plot some summaries from the DE analysis
    """
    input:
        DEGs_sc=f"results/DE_analysis/DEGs_sc.csv",
        DEGs_pb=f"results/DE_analysis/DEGs_pb.csv",
        DEGs_both=f"results/DE_analysis/DEGs_both.csv",
    output:
        overlap=f"results/DE_analysis/plots/sc_pb_overlap.pdf",
        FDR_correlation=f"results/DE_analysis/plots/FDR_correlation.pdf",
        log2FC_correlation=f"results/DE_analysis/plots/log2FC_correlation.pdf",
    log:
        "logs/DE_analysis/DEA_plots.log",
    conda:
        "../envs/DE_analysis_plot.yaml"
    threads: config["DE_analysis"]["DEA_plots"]["cores"]
    resources:
        mem_mb=config["DE_analysis"]["DEA_plots"]["mem_mb"],
        time=config["DE_analysis"]["DEA_plots"]["time"],
        queue=config["queues"]["cpu"],
    message:
        "Plot summary of DE analysis."
    script:
        "../scripts/DE_analysis/DEA_plots.R"
