configfile: "config/config.yaml"


rule annotate_and_save:
    """
    Modify the metadata layer of the scRNAseq dataset to include the cell type annotation
    from cytopus, as specified on the user-provided cell type conversion dictionary.
    Save the dataset in h5ad, 10XGenomics mtx format and rds.
    """
    input:
        sc_dataset = config["scRNAseq"],
        dictionary = config["preprocess"]["annotate_and_save"]["celltype_conversion_dictionary"]
    output:
        anndata  = "results/preprocess/AnnData/"+config["analysis_name"]+".h5ad",
        rds      = "results/preprocess/rds/"+config["analysis_name"]+".rds",
        matrix   = "results/preprocess/mtx/"+config["analysis_name"]+"_matrix.mtx",
        barcodes = "results/preprocess/mtx/"+config["analysis_name"]+"_barcodes.tsv",
        genes    = "results/preprocess/mtx/"+config["analysis_name"]+"_genes.tsv"
    log:
        "logs/preprocess/annotate_and_save.log"
    message:
        "Annotating Seurat object's metadata with cytopus cell types and saving as .rds, 10x-Genomics-formatted mtx and h5ad ..."
    threads: 1
    resources: 
        mem_mb = config["preprocess"]["annotate_and_save"]["rstudio_memory"]+500,
        time   = config["preprocess"]["annotate_and_save"]["time"]
    script:
        "scripts/preprocess/annotate_and_save.R"