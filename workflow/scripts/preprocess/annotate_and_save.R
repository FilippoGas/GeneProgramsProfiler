# ==============================================================================
# Script: annotate_and_save.R
# Author: Filippo Gastaldello
# Date: 16/07/2026
# Description: 
#   Match the cell type annotation present in the scRNAseq dataset with the cell
#   types available in cytopus. Save the dataset in the formats needed in the 
#   following steps of the pipeline (rds, 10XGenomics  mtx format and h5ad)
#
# Snakemake Expected Inputs:
#   - snakemake@input[["sc_dataset"]] : Path to seurat object (.rds)
#   - snakemake@input[["dictionary"]] : Path to cell type conversion dictionary
#                                       (original annotation to cytopus)
#
# Snakemake Expected Outputs:
#   - snakemake@output[["anndata"]]  : Path to save expression matrix (.h5ad)
#   - snakemake@output[["rds"]]      : Path to save seurat object (.rds)
#   - snakemake@output[["matrix"]]   : Path to save expression matrix (.mtx)
#   - snakemake@output[["barcodes"]] : Path to save cell IDs (.mtx)
#   - snakemake@output[["genes"]]    : Path to save gene names (.mtx)
#
# Snakemake Expected Params:
#   - snakemake@params[["annotation_column"]] : Name of the column containing
#                                               cell type annotation in the
#                                               original dataset
# ==============================================================================

# Setup Logging ----------------------------------------------------------------
# Redirect all output and messages to the Snakemake log file
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# Load Libraries ---------------------------------------------------------------
suppressPackageStartupMessages({
        library(tidyverse)
        library(jsonlite)
        library(Matrix)
        library(Seurat)
        library(SingleCellExperiment)
        library(zellkonverter)
})

message("Starting R script \"annotate_and_save.R\"...")

# 1. Load Data ####
# ------------------------------------------------------------------------------
message("Loading Seurat object: ", snakemake@input[["sc_dataset"]])
data <- read_rds(snakemake@input[["sc_dataset"]])
message("Done")
message("Loading cell type conversion dictionary: ",
        snakemake@input[["dictionary"]])
celltype_conversion_dict <- read_json(snakemake@input[["dictionary"]],
                                      show_col_types = FALSE)
celltype_conversion_dict <- data.frame(
                        "celltype"=names(celltype_conversion_dict),
                        "cytopus"=unname(unlist(celltype_conversion_dict))
                                       )
message("Done")

# 2. Load Parameters ####
# ------------------------------------------------------------------------------
message("Loading annotation column name: ",
        snakemake@params[["annotation_colname"]])
celltype_annotation_colname <- snakemake@params[["annotation_colname"]]
message("Done")

# 3. Add cytopus cell type annotation ####
# ------------------------------------------------------------------------------
message("Adding cell type annotations from Cytopus... ")
# Make sure the default assay for the seurat object is set to "RNA"
DefaultAssay(data) <- "RNA"
# Rename original celltype annotation column to "celltype"
data@meta.data <- data@meta.data %>%
                        dplyr::rename("celltype" = celltype_annotation_colname)
# Add annotations to metadata in all datasets
data@meta.data <- data@meta.data %>% 
        rownames_to_column("barcode") %>%
        left_join(cytopus_dict, by = "celltype") %>%
        column_to_rownames("barcode")
message("Done")

# 4. Save Outputs ####
# ------------------------------------------------------------------------------
message("Saving dataset as h5ad to: ", snakemake@output[["anndata"]])
writeH5AD(
        as.SingleCellExperiment(data),
        file = snakemake@output[["anndata"]]
)
message("Done")
message("Saving dataset as rds to ", snakemake@output[["rds"]])
saveRDS(data, snakemake@output[["rds"]])
message("Done")
message("Saving count matrix as .mtx to: ", snakemake@output[["matrix"]])
counts <- data@assays$RNA$counts
writeMM(
        counts,
        snakemake@output[["matrix"]]
)
message("Done")
message("Saving cell barcodes as .mtx to: ", snakemake@output[["barcodes"]])
write.table(
        as.data.frame(colnames(counts)),
        snakemake@output[["barcodes"]],
        col.names = FALSE,
        row.names = FALSE,
        sep = "\t"
)
message("Done")
message("Saving gene names as .mtx to: ", snakemake@output[["genes"]])
features <- data.frame(
                "gene_id"    = rownames(count),
                "gene_names" = rownames(count),
                type = "Gene Expression"
)
write.table(
        as.data.frame(features),
        snakemake@output[["genes"]]
)
message("Done")
# Close Logging ----------------------------------------------------------------
sink(type="message")
sink()