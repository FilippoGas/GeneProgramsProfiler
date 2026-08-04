# ==============================================================================
# Script: run_ora.R
# Author: Filippo Gastaldello
# Date: 04/08/2026
# Description: 
#   Performs Over Representation Analysis on the results of the differential
#   analyses
#
# Snakemake Expected Inputs:
#   - snakemake@input[["DEGs_both"]] : Path to results of DE analyses (from the
#                                      merged csv both sc and pb can be 
#                                      retrieved)
#   - snakemake@input[["cytopus]] : Path to the json dictionary with all cytopus
#                                   gene sets used for the enrichment
#   - snakemake@input[["gene_list]] : List of all genes present in the dataset,
#                                     it will be used as background during ORA
#
# Snakemake Expected Outputs:
#   - snakemake@output[[1]] : path to fake output used to signal that the rule
#                             has finished. Real outputs will be saved in the 
#                             same folder with name generated dynamically.
#
# Snakemake Expected Params:
#   - snakemake@params[["padj_thresh"]] : Value to use as threshold for the 
#                                         adjusted pvalue in the ORA results
# ==============================================================================

# Setup Logging ----------------------------------------------------------------
# Redirect all output and messages to the Snakemake log file
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# Load Libraries ---------------------------------------------------------------
suppressPackageStartupMessages({
        library(tidyverse)
        library(ggplot2)
        library(jsonlite)
        library(clusterProfiler)
        library(AnnotationDbi)
        library(org.Hs.eg.db)
})
# Prevent R from generating the default Rplots.pdf file
pdf(NULL)

message("Starting R script \"run_ora.R\"...")

# 1. Load Data ####
# ------------------------------------------------------------------------------
message("Loading results from differential analysis from ",
        snakemake@input[["DEGs_both"]])
DEGs_both <- read_csv(snakemake@input[["DEGs_both"]])
message("Loading cytopus gene sets from ", snakemake@input[["cytopus"]])
cytopus <- read_json(snakemake@input[["cytopus"]])
message("Done")
message("Loading gene list from ",snakemake@input[["gene_list"]])
gene_list <- readLines(snakemake@input[["gene_list"]])
message("Done")
message("Loading padj threshold ",snakemake@params[["padj_thresh"]])
padj_thresh <- snakemake@params[["padj_thresh"]]
message("Done")

# 2. Prepare cytopus gene sets to be used in ClusterProfiler ####
# ------------------------------------------------------------------------------
message("Preparing cytopus gene sets to be used in Cluster Profiler ...")
analyses         <- c("sc", "pb", "both")
class_columns    <- list("sc"="class_sc",
                         "pb"="class_pb",
                         "both"="class_pb")
sig_columns      <- list("sc"="is_sig_sc",
                         "pb"="is_sig_pb",
                         "both"="both_sig")
# Get all valid human gene symbols
valid_symbols <- AnnotationDbi::keys(org.Hs.eg.db, keytype = "SYMBOL")
# Prepare background universe set (all genes that were tested in DE analysis)
universe_symbol <- intersect(gene_list, valid_symbols)
message("Done")

# 3. Run ORA for single cell, pseudobulk and common DEGs ####
# ------------------------------------------------------------------------------
message("Running ORA ...")
# Read output path
base_out_dir <- str_split_i(snakemake@output[[1]], ".done", 1)
for (cell_type in unique(DEGs_both$celltype)) {
        # Prepare Cytopus term2gene once per cell type
        cytopus_list <- c(cytopus[[cell_type]], cytopus[["global"]])
        cytopus_list <- lapply(cytopus_list, function(x) as.character(unlist(x, use.names = FALSE)))
        cytopus_term2gene <- data.frame(
                term = rep(names(cytopus_list), lengths(cytopus_list)),
                gene = unlist(cytopus_list, use.names = FALSE)
        )
        for (analysis in analyses) {
                # Loop through both UP and DOWN directions
                for (dir in c("UP", "DOWN")) {
                        # Extract significant DEGs for the current analysis type and direction
                        sig_genes <- DEGs_both %>% 
                                dplyr::filter(.data[[sig_columns[[analysis]]]] == TRUE,
                                              .data[[class_columns[[analysis]]]] == dir,
                                              celltype == cell_type) %>% 
                                dplyr::pull(gene)
                        # Keep only the genes that actually exist in the database
                        mapped_genes <- intersect(sig_genes, valid_symbols)
                        # Skip to the next iteration if there are no significant genes for this set
                        if (length(mapped_genes) == 0) next 
                        # Run ORA against Cytopus
                        ora_results <- enricher(mapped_genes,
                                                minGSSize = 5,
                                                universe  = universe_symbol,
                                                TERM2GENE = cytopus_term2gene)
                        # Check if the result object is valid and contains significant adjusted p-values
                        if (!is.null(ora_results) && nrow(ora_results@result) > 0 && min(ora_results@result$p.adjust, na.rm = TRUE) < padj_thresh) {
                                # Format the title text based on direction
                                dir_text <- ifelse(dir == "UP", "up-regulated", "down-regulated")
                                p <- dotplot(ora_results, showCategory=20, font.size=10, label_format=70) +
                                        scale_size_continuous(range=c(1, 7)) +
                                        theme_minimal() +
                                        ggtitle(paste0("Cytopus Enrichment of ", dir_text, " genes for ", cell_type, " cells."))
                                
                                # Define output directory and filename prefix
                                out_dir     <- paste0(base_out_dir, analysis, "/")
                                file_prefix <- paste0(out_dir, gsub("/", "_", cell_type), "_", dir)
                                
                                ggsave(p, filename = paste0(file_prefix, ".pdf"), device = "pdf", width = 12, height = 8, create.dir = TRUE)
                                write_csv(ora_results@result, file = paste0(file_prefix, ".csv"))
                        }
                }
        }
}
write("Done.", file = snakemake@output[[1]])
message("Done")

# Close Logging ----------------------------------------------------------------
sink(type="message")
sink()