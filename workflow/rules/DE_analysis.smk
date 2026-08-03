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
        DEGs_pb=f"results/DE_analysis/DEGs_pb",
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
