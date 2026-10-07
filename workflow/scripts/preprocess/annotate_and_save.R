# ==============================================================================
# Script: annotate_and_save.R
# Author: Filippo Gastaldello
# Date: 16/07/2026
# Description: 
#   Match the cell type annotation present in the scRNAseq dataset with the cell
#   types available in cytopus. Save the dataset in the formats needed in the 
#   following steps of the pipeline (rds, 10XGenomics  mtx format and h5ad)
#
#   Memory strategy: when the count/data layers are BPCells IterableMatrix
#   objects (on-disk), every export streams through the matrix in bounded-size
#   blocks and NEVER materializes the full matrix in memory:
#     - rds  : BPCells directories are copied (chunked) into the results dir
#              and the Seurat layers are re-pointed at the copies
#     - mtx  : a column-block writer replaces Matrix::writeMM (which has no
#              method for IterableMatrix) and emits standard MatrixMarket
#     - h5ad : the X layers are written with BPCells write_matrix_anndata_hdf5
#              (chunked); obs/var/root metadata are attached via anndataR's
#              rhdf5-based writers (no python/anndata needed)
#   In-memory (dgCMatrix) inputs keep the original byte-compatible behaviour.
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
        library(BPCells)
        library(rhdf5)
})

# Helper functions -------------------------------------------------------------
# Stream a BPCells (or in-memory) matrix to MatrixMarket (.mtx) in column blocks
# so the full matrix is never materialized in memory.
is_iterable_on_disk <- function(mat) {
        inherits(mat, "IterableMatrix")
}

write_mtx_streamed <- function(mat, path, block_cols = 4096L) {
        nrow <- nrow(mat)
        ncol <- ncol(mat)
        # dgCMatrix block iteration is column-major, so .mtx is written in
        # column-major order like Matrix::writeMM(., dgCMatrix).
        block_seq <- seq(1L, ncol, by = block_cols)
        # Determine MatrixMarket value type: "integer" if every stored value is
        # a whole number, "real" otherwise (mirrors Matrix::writeMM).
        field <- "integer"
        con <- file(path, open = "wt")
        on.exit(close(con), add = TRUE)
        # Need nnz for the header line before streaming values -> two passes.
        nnz <- 0
        for (start in block_seq) {
                end <- min(start + block_cols - 1L, ncol)
                blk <- as(mat[, start:end], "dgCMatrix")
                nnz <- nnz + length(blk@x)
                if (field == "integer" && any(blk@x != floor(blk@x))) {
                        field <- "real"
                }
        }
        writeLines(
                c(
                        sprintf("%s", paste0(
                                "%%MatrixMarket matrix coordinate ", field,
                                " general"
                        )),
                        paste(
                                format(c(nrow, ncol, nnz), scientific = FALSE,
                                       trim = TRUE),
                                collapse = " "
                        )
                ),
                con
        )
        for (start in block_seq) {
                end <- min(start + block_cols - 1L, ncol)
                blk <- as(mat[, start:end], "dgCMatrix")
                if (length(blk@x) == 0) {
                        next
                }
                # dgC stores entries column-major: @i gives the (0-based) row of
                # each entry; the column is recovered from the column pointers.
                i <- blk@i + 1L
                j <- rep.int(
                        seq.int(start, end),
                        diff(blk@p)
                )
                if (field == "integer") {
                        vals <- format(blk@x, scientific = FALSE, trim = TRUE)
                } else {
                        vals <- format(blk@x, scientific = FALSE, digits = 16,
                                       trim = TRUE)
                }
                lines <- paste(i, j, vals, sep = " ")
                writeLines(lines, con)
        }
        invisible(path)
}

# Copy an on-disk BPCells layer into the results tree (chunked) and return an
# IterableMatrix pointing at the copy (used to make the saved rds self-contained).
copy_bpcells_layer <- function(mat, dest_dir) {
        dir.create(dirname(dest_dir), recursive = TRUE, showWarnings = FALSE)
        if (dir.exists(dest_dir)) {
                unlink(dest_dir, recursive = TRUE)
        }
        write_matrix_dir(mat, dest_dir)
        open_matrix_dir(dest_dir)
}

# Write an .h5ad whose expression matrices are streamed by BPCells and whose
# obs/var/root metadata is (re)written with anndataR's rhdf5-based writers.
# mat_layers: named list of IterableMatrix, names become layer names.
write_h5ad_streamed <- function(mat_layers, path, obs_df, var_df) {
        if (file.exists(path)) {
                unlink(path)
        }
        # 1. Stream all expression layers into the file (BPCells writes them as
        #    anndata-compatible csr_matrix groups, transposed to cells x genes).
        for (layer_name in names(mat_layers)) {
                # The returned object holds an open handle on the file; discard
                # it explicitly so rhdf5 can open the file for writing below.
                res <- write_matrix_anndata_hdf5(
                        mat_layers[[layer_name]],
                        path,
                        group = paste0("layers/", layer_name)
                )
                rm(res)
        }
        gc(verbose = FALSE)
        # 2. Replace BPCells' minimal /obs and /var (dimnames only) with the
        #    full dataframes, and add the root anndata encoding attributes.
        h5 <- rhdf5::H5Fopen(path, flags = "H5F_ACC_RDWR")
        on.exit(rhdf5::H5Fclose(h5), add = TRUE)
        for (grp in c("obs", "var")) {
                tryCatch(
                        rhdf5::h5delete(h5, paste0("/", grp)),
                        error = function(e) NULL
                )
        }
        anndataR:::hdf5_write_attribute(h5, "/", "encoding-type", "anndata")
        anndataR:::hdf5_write_attribute(h5, "/", "encoding-version", "0.1.0")
        anndataR:::write_h5ad_element(obs_df, h5, "obs")
        anndataR:::write_h5ad_element(var_df, h5, "var")
        invisible(path)
}

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
raw_counts <- data@assays$RNA$counts
on_disk <- is_iterable_on_disk(raw_counts)
if (on_disk) {
        message("Detected on-disk (BPCells) integrity layers: streaming all exports.")
}

# 1. Save RDS 
# ------------------------------------------------------------------------------
message("Saving dataset as rds to ", snakemake@output[["rds"]])
if (on_disk) {
        # Copy the BPCells layers (chunked) into the results tree and re-point
        # the object at the copies so the saved rds is fully self-contained.
        rds_bp_dir <- file.path(dirname(snakemake@output[["rds"]]), "bpcells")
        counts_layer_old <- data@assays$RNA$counts
        counts_layer_new <- copy_bpcells_layer(
                counts_layer_old,
                file.path(rds_bp_dir, "counts")
        )
        data <- SetAssayData(
                data,
                assay = "RNA",
                layer = "counts",
                new.data = counts_layer_new
        )
        if ("data" %in% Layers(data)) {
                data_layer_old <- data@assays$RNA$data
                if (is_iterable_on_disk(data_layer_old)) {
                        data_layer_new <- copy_bpcells_layer(
                                data_layer_old,
                                file.path(rds_bp_dir, "data")
                        )
                        data <- SetAssayData(
                                data,
                                assay = "RNA",
                                layer = "data",
                                new.data = data_layer_new
                        )
                }
        }
}
saveRDS(data, snakemake@output[["rds"]])
message("Done")

# 2. Save MTX and features using the original raw counts
# ------------------------------------------------------------------------------
message("Saving raw count matrix as .mtx to: ", snakemake@output[["matrix"]])
if (on_disk) {
        write_mtx_streamed(raw_counts, snakemake@output[["matrix"]])
} else {
        Matrix::writeMM(
                raw_counts,
                snakemake@output[["matrix"]]
        )
}
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
# ------------------------------------------------------------------------------
message("Preparing matrix for h5ad export...")
# Find variable features and subset before h5ad conversion
# FIXME: make nfeatures configurable via config.yaml
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 4000) 
# Subset the data to only keep those genes
data <- data[VariableFeatures(data), ]
log_norm_matrix <- GetAssayData(data, assay = "RNA", layer = "data")

message("Saving log-normalized dataset as h5ad to: ",
        snakemake@output[["anndata"]])
if (on_disk) {
        # The in-memory behaviour below stores BOTH the "counts" and "data"
        # layers as the log-normalized variable-feature matrix; reproduce that
        # exactly, but stream the values instead of densifying them.
        write_h5ad_streamed(
                mat_layers = list(
                        "counts" = log_norm_matrix,
                        "data" = log_norm_matrix
                ),
                path = snakemake@output[["anndata"]],
                obs_df = as.data.frame(data@meta.data),
                var_df = data[["RNA"]][[]]
        )
} else {
        data <- SetAssayData(
                data,
                assay = "RNA",
                layer = "counts",
                new.data = log_norm_matrix)
        write_h5ad(
                data,
                path = snakemake@output[["anndata"]]
        )
}
message("Done")
# Close Logging ----------------------------------------------------------------
sink(type="message")
sink()