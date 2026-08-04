rule prepare_cytopus_list:
    """
    Download genesets from cytopus for the cell types present in the scRNAseq dataset.
    For more details visit cytopus (https://github.com/wallet-maker/cytopus)
    """
    input:
        config["celltype_conversion_dictionary"],
    output:
        f"results/{config["analysis_name"]}/spectra/cytopus_Gene_sets.json",  #TODO move analysis name after results/ and do the same for all rules. so results will contain all differen analysis
    log:
        f"logs/{config["analysis_name"]}/spectra/prepare_cytopus_list.log",
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
