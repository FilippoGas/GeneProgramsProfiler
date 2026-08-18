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
#   - snakemake@params[["sample_col"]] : Name of the column holding sample 
#                                        names in the metadata
#   - snakemake@params[["condition_col"]] : Name of the column holding sample 
#                                        condition in the metadata
#   - snakemake@params[["phase_col"]] : Name of the column holding cell cycle 
#                                       phase in the metadata
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
        library(anndataR)
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
message("Loading diagnosis column name: ",
        snakemake@params[["condition_col"]])
condition_annotation_colname <- snakemake@params[["condition_col"]]
message("Done")
message("Loading sample name column name: ",
        snakemake@params[["sample_col"]])
sample_name_colname <- snakemake@params[["sample_col"]]
message("Done")
message("Loading cell cycle phase column name: ",
        snakemake@params[["phase_col"]])
cell_cycle_phase_colname <- snakemake@params[["phase_col"]]
message("Done")

# 3. Add cytopus cell type annotation ####
# ------------------------------------------------------------------------------
message("Adding cell type annotations from Cytopus... ")
# Make sure the default assay for the seurat object is set to "RNA"
DefaultAssay(data) <- "RNA"
message(colnames(data@meta.data))
# Rename original metadata columns
data@meta.data <- data@meta.data %>%
                        dplyr::rename(
                                "celltype" = all_of(
                                        celltype_annotation_colname
                                        )
                                )
data@meta.data <- data@meta.data %>%
                        dplyr::rename(
                                "diagnosis" = all_of(
                                        condition_annotation_colname
                                )
                        )
data@meta.data <- data@meta.data %>%
        dplyr::rename(
                "sample_name" = all_of(
                        sample_name_colname
                )
        )
data@meta.data <- data@meta.data %>%
        dplyr::rename(
                "phase" = all_of(
                        cell_cycle_phase_colname
                )
        )
# Add annotations to metadata 
data@meta.data <- data@meta.data %>% 
        rownames_to_column("barcode") %>%
        left_join(celltype_conversion_dict, by = "celltype") %>%
        column_to_rownames("barcode")
message("Done")

# 4. Save Outputs ####
# ------------------------------------------------------------------------------
# 1. Save RDS 
message("Saving dataset as rds to ", snakemake@output[["rds"]])
saveRDS(data, snakemake@output[["rds"]])
message("Done")

# 2. Save MTX and features using the original raw counts
message("Saving raw count matrix as .mtx to: ", snakemake@output[["matrix"]])
raw_counts <- data@assays$RNA$counts
writeMM(
        raw_counts,
        snakemake@output[["matrix"]]
)
message("Done")

message("Saving cell barcodes as .mtx to: ", snakemake@output[["barcodes"]])
write.table(
        as.data.frame(colnames(raw_counts)),
        snakemake@output[["barcodes"]],
        col.names = FALSE,
        row.names = FALSE,
        quote = FALSE,
        sep = "\t"
)
message("Done")

message("Saving gene names as .mtx to: ", snakemake@output[["genes"]])
features <- data.frame(
        "gene_id"    = rownames(raw_counts),
        "gene_names" = rownames(raw_counts),
        "type"       = "Gene Expression"
)
write.table(
        features,
        snakemake@output[["genes"]],
        sep = "\t",
        row.names = FALSE,
        col.names = FALSE,
        quote = FALSE
)
message("Done")

# 3. Modify the object strictly for the h5ad converter
message("Preparing matrix for h5ad export...")
# Find variable features and subset before h5ad conversion
#TODO make this optional and configurable from the config file
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 4000) 
# Subset the data to only keep those genes
data <- data[VariableFeatures(data), ]
log_norm_matrix <- GetAssayData(data, assay = "RNA", layer = "data")
data <- SetAssayData(
                data,
                assay = "RNA",
                layer = "counts",
                new.data = log_norm_matrix)

message("Saving log-normalized dataset as h5ad to: ",
        snakemake@output[["anndata"]])
write_h5ad(
        data,
        path = snakemake@output[["anndata"]]
)
message("Done")
# Close Logging ----------------------------------------------------------------
sink(type="message")
sink()