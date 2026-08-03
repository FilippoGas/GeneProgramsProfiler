# ==============================================================================
# Script: run_DE_analysis.R
# Author: Filippo Gastaldello
# Date: 03/08/2026
# Description: 
#   Run differential expression analysis on the input dataset.
#   Performs analysis in single cell and pseudo bulk modes.
#
# Snakemake Expected Inputs:
#   - snakemake@input[["sc_dataset"]] : Path to seurat object (.rds)
#
# Snakemake Expected Outputs:
#   - snakemake@output[["DEGs_sc"]]   : Path to save DEGs according to single 
#                                       cell analysis
#   - snakemake@output[["DEGs_pb"]]   : Path to save DEGs according to 
#                                       pseudobulk analysis
#   - snakemake@output[["DEGs_both"]] : Path to save DEGs according to both 
#                                       analysis
# Snakemake expected params:
#   - snakemake@params[["case"]]   : Name of case condition
#   - snakemake@params[["control]] : Name of control condition
# ==============================================================================

# Setup Logging ----------------------------------------------------------------
# Redirect all output and messages to the Snakemake log file
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# Load Libraries ---------------------------------------------------------------
suppressPackageStartupMessages({
        library(tidyverse)
        library(Seurat)
        library(parallel)
})

message("Starting R script \"annotate_and_save.R\"...")

# 1. Load Data ####
# ------------------------------------------------------------------------------
message("Loading Seurat object: ", snakemake@input[["sc_dataset"]])
data <- read_rds(snakemake@input[["sc_dataset"]])
message("Done")
message("Loading params ...")
message("Case condition: ", snakemake@params[["case"]])
case <- snakemake@params[["case"]]
message("Control condition: ", snakemake@params[["control"]])
control <- snakemake@params[["control"]]
message("FDR threshold: ", snakemake@params[["FDR_thresh"]])
fdr_thresh <- snakemake@params[["FDR_thresh"]]
message("logFC thresh: ", snakemake@params[["logFC_thresh"]])
logFC_thresh <- snakemake@params[["logFC_thresh"]]
message("Sample name column: ", snakemake@params[["sample_col"]])
sample_col <- snakemake@params[["sample_col"]]
message("Condition name column: ", snakemake@params[["condition_col"]])
condition_col <- snakemake@params[["condition_col"]]
message("Done")

# 2. Run single cell differential analysis ####
# ------------------------------------------------------------------------------
message("Running single cell differential analysis ...")
# Create new metadata column to use as comparison
data@meta.data <- data@meta.data %>% 
        mutate(celltype_condition = paste0(
                                        celltype,
                                        "-",
                                        Diagnosis)
               )
Idents(data) <- "celltype_condition"
# Run DE for each celltype
res_DE <- mclapply(unique(data@meta.data[,"celltype"]),
                   function(cell_type){
                           return(
                                FindMarkers(data,
                                            ident.1 = paste0(cell_type,
                                                             "-",
                                                             case),
                                            ident.2 = paste0(cell_type,
                                                             "-",
                                                             control)
                                              ) %>%
                                mutate(celltype = cell_type,
                                       FDR = p.adjust(p_val,
                                                      method = "fdr")) %>%
                                rownames_to_column(var = "gene")
                                )
                   },
                   mc.cores = snakemake@threads-1,
                   mc.cleanup = TRUE)
# Combine results from different cell types
DEGs_sc <- bind_rows(res_DE)
# Retrieve smallest non-zero pval
lowest_p <- (DEGs_sc$FDR %>% unique() %>% sort())[2]
# Label deregulated genes and add rank metric (-log10(FDR)*log2FC)
# lowest_p in used as a pseudo-count during the rank computation to avoid
# log(0) (Some pvaules might be zero)
DEGs_sc <- DEGs_sc %>% 
        mutate(is_sig = ifelse(FDR < fdr_thresh & abs(avg_log2FC) >logFC_thresh,
                               TRUE,
                               FALSE),
               class  = ifelse(avg_log2FC > 0,
                               "UP",
                               "DOWN"),
               rank   = -log10(FDR+lowest_p)*avg_log2FC)
message("Done")

# 3. Run pseudobulk differential analysis ####
# ------------------------------------------------------------------------------
message("Running pseudo-bulk differential analysis ...")
# Aggregate expression per sample
data_pb <- AggregateExpression(data,
                               assays = "RNA",
                               return.seurat = TRUE,
                               group.by = c(sample_col,
                                            condition_col,
                                            "celltype")
                               )
# Create new idents
data_pb$celltype_condition <- paste0(data_pb@meta.data[, "celltype"],
                                     "-",
                                     data_pb@meta.data[, condition_col])
Idents(data_pb) <- "celltype_condition"
# Run DE for each cell type
res_DE_pb <- lapply(unique(data_pb@meta.data[, "celltype"]),
                    function(cell_type){
                            return(FindMarkers(data_pb,
                                               ident.1 = paste0(cell_type,
                                                                "-",
                                                                case),
                                               ident.2 = paste0(cell_type,
                                                                "-",
                                                                control),
                                               test.use = "DESeq2") %>%
                                   mutate(celltype = cell_type,
                                          FDR = p.adjust(p_val,
                                                         method = "fdr")) %>%
                                   rownames_to_column(var = "gene")
                            )
                    }
)
# Combine results from different cell types
DEGs_pb <- bind_rows(res_DE_pb)
# Retrieve smallest non-zero pval
lowest_p <- (DEGs_pb$FDR %>% unique() %>% sort())[2]
# Label deregulated genes and add rank metric (-log10(FDR)*log2FC)
# lowest_p in used as a pseudo-count during the rank computation to avoid
# log(0) (Some pvaules might be zero)
DEGs_pb <- DEGs_pb %>%
        mutate(is_sig = ifelse(FDR < fdr_thresh & abs(avg_log2FC) >logFC_thresh,
                               TRUE,
                               FALSE),
               is_sig = ifelse(is.na(p_val),
                               FALSE,
                               is_sig),
               class  = ifelse(avg_log2FC > 0,
                               "UP",
                               "DOWN"),
               rank   = -log10(FDR+lowest_p)*avg_log2FC)
message("Done")

# 4. Save DEGs csv tables ####
# ------------------------------------------------------------------------------
message("Saving single cell DEGs to ", snakemake@output[["DEGs_sc"]], " ...")
write_csv(DEGs_sc, file = snakemake@output[["DEGs_sc"]])
message("Done")
message("Saving pseudo-bulk DEGs to ", snakemake@output[["DEGs_pb"]], " ...")
write_csv(DEGs_pb, file = snakemake@output[["DEGs_pb"]])
message("Done")
message("Merging single cell and pseudo-bulk DEGs ...")
DEGs_both <- DEGs_sc %>% full_join(DEGs_pb,
                                    by = c("gene", "celltype"),
                                    suffix = c("_sc", "_pb"))
# Create flag for singificance both in sc and pb
DEGs_both <- DEGs_both %>%
        mutate(both_sig = ifelse(is_sig_sc == TRUE & is_sig_pb == TRUE,
                                 TRUE,
                                 FALSE)
               )
message("Done")
message("Saving merged DEGs to ", snakemake@output[["DEGs_both"]], " ...")
write_csv(DEGs_both, file = snakemake@output[["DEGs_both"]])
message("Done")

# Close Logging ----------------------------------------------------------------
sink(type="message")
sink()