rule annotate_and_save:
    """
    Modify the metadata layer of the scRNAseq dataset to include the cell type annotation
    from cytopus, as specified on the user-provided cell type conversion dictionary.
    Save the dataset in h5ad, 10XGenomics mtx format and rds.
    """
    input:
        sc_dataset=config["scRNAseq"],
        dictionary=config["celltype_conversion_dictionary"],
    output:
        anndata=f"results/{config["analysis_name"]}/preprocess/AnnData/adata.h5ad",
        rds=f"results/{config["analysis_name"]}/preprocess/rds/seurat.rds",
        matrix=f"results/{config["analysis_name"]}/preprocess/mtx/matrix.mtx",
        barcodes=f"results/{config["analysis_name"]}/preprocess/mtx/barcodes.tsv",
        genes=f"results/{config["analysis_name"]}/preprocess/mtx/genes.tsv",
    log:
        f"logs/{config["analysis_name"]}/preprocess/annotate_and_save.log",
    conda:
        "../envs/preprocess.yaml"
    threads: config["preprocess"]["annotate_and_save"]["cores"]
    resources:
        mem_mb=config["preprocess"]["annotate_and_save"]["rstudio_memory"] + 500,
        time=config["preprocess"]["annotate_and_save"]["time"],
        queue=config["queues"]["cpu"],
    params:
        annotation_colname=config["preprocess"]["annotate_and_save"][
            "celltype_annotation_colname"
        ],
        sample_col=config["preprocess"]["annotate_and_save"]["sample_column"],
        condition_col=config["preprocess"]["annotate_and_save"]["condition_column"],
        phase_col=config["preprocess"]["annotate_and_save"]["cell_cycle_phase_column"],
    message:
        "Annotate Seurat object's metadata with cytopus cell types and save as .rds, 10x-Genomics-formatted mtx and h5ad."
    script:
        "../scripts/preprocess/annotate_and_save.R"
