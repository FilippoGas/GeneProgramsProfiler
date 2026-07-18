configfile: ".test/config/config.yaml"


rule prepare_cytopus_list:
    """
    Download genesets from cytopus for the cell types present in the scRNAseq dataset.
    For more details visit cytopus (https://github.com/wallet-maker/cytopus)
    """
    input:
        config["celltype_conversion_dictionary"],
    output:
        f"results/spectra/{config["analysis_name"]}_cytopus_Gene_sets.json"
    log:
        "logs/spectra/prepare_cytopus_list.log",
    conda:
        "../envs/spectra.yaml"
    threads: 
        config["spectra"]["prepare_cytopus_list"]["cores"]
    resources:
        mem_mb=config["spectra"]["prepare_cytopus_list"]["mem_mb"],
        time=config["spectra"]["prepare_cytopus_list"]["time"],
        queue=config["queues"]["cpu"]
    params:
        global_celltype=config["spectra"]["prepare_cytopus_list"]["global_celltype"]
    message:
        "Retrieve cytopus gene sets for the cell types in the dataset."
    script:
        "../scripts/spectra/prepare_cytopus_list.py"

rule run_spectra:
    """
    Run spectra to quantify gene programs activation in the scRNAseq dataset
    """
    input:
    output:
    log:
    conda:
    threads