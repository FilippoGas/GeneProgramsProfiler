# ==============================================================================
# Script: run_gsea.R
# Author: Filippo Gastaldello
# Date: 04/08/2026
# Description: 
#   Performs Gene Set Enrichment Analysis on the results of the differential
#   analyses
#
# Snakemake Expected Inputs:
#   - snakemake@input[["DEGs_both"]] : Path to results of DE analyses (from the
#                                      merged csv both sc and pb can be 
#                                      retrieved)
#   - snakemake@input[["cytopus]] : Path to the json dictionary with all cytopus
#                                   gene sets used for the enrichment
#
# Snakemake Expected Outputs:
#   - snakemake@output[[1]] : path to fake output used to signal that the rule
#                             has finished. Real outputs will be saved in the 
#                             same folder with name generated dynamically.
#
# Snakemake Expected params:
#   - snakemake@params[["padj_thresh"]] : Value to use as threshold for the 
#                                         adjusted pvalue in the fgsea results
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
        library(ggplotify)
        library(patchwork)
        library(msigdbr)
        library(jsonlite)
        library(fgsea)
})

message("Starting R script \"run_gsea.R\"...")

# 1. Load Data ####
# ------------------------------------------------------------------------------
message("Loading results from differential analysis from ",
        snakemake@input[["DEGs_both"]])
DEGs_both <- read_csv(snakemake@input[["DEGs_both"]])
message("Loading cytopus gene sets from ",
        snakemake@input[["cytopus"]])
cytopus <- read_json(snakemake@input[["cytopus"]])
message("Done")
message("Loading padj threshold: ",
        snakemake@params[["padj_thresh"]])
padj_thresh <- snakemake@params[["padj_thresh"]]
message("Done")

# 2. Prepare cytopus gene sets for enrichment with fgsea ####
# ------------------------------------------------------------------------------
message("Reformatting cytopus gene sets to be used in fgsea ...")
cytopus_list        <- unlist(cytopus, recursive = FALSE)
names(cytopus_list) <- str_split_i(names(cytopus_list), "\\.",2) 
cytopus_list        <- cytopus_list[!duplicated(names(cytopus_list))]
cytopus_list        <- lapply(cytopus_list, function(x){unlist(x)})
message("Done")

# 3. Run GSEA for single cell and pseudobulk DEGs ####
# ------------------------------------------------------------------------------
message("Running GSEA ...")
# Read output path
out_dir <- str_split_i(snakemake@output[[1]], ".done", 1)
# Define configurations for each analysis
analysis_configs <- list(
        sc = list(
                sig_col    = "is_sig_sc",
                rank_col   = "rank_sc",
                title_desc = "single cell DE analysis",
                out_folder = "sc"
        ),
        pb = list(
                sig_col    = "is_sig_pb",
                rank_col   = "rank_pb",
                title_desc = "pseudo-bulk DE analysis",
                out_folder = "pb"
        )
)
# For each celltype run gsea for sc and pb DEGs 
for (cell_type in unique(DEGs_both$celltype)) {
        # Iterate over sc and pb analyses
        for (analysis_name in names(analysis_configs)) {
                
                # Get configuration for the current analysis
                cfg <- analysis_configs[[analysis_name]]
                
                # Create rank for this celltype based on the configuration
                celltype_DEGs <- DEGs_both %>% 
                        dplyr::filter(celltype == cell_type,
                                      !is.na(.data[[cfg$sig_col]])) %>% 
                        dplyr::arrange(desc(.data[[cfg$rank_col]]))
                
                rank        <- celltype_DEGs[[cfg$rank_col]]
                names(rank) <- celltype_DEGs$gene
                rank        <- rank[!is.na(rank)]
                set.seed(42)
                fgsea_results <- fgsea(pathways = cytopus_list, 
                                       stats    = rank,
                                       minSize  = 5,
                                       maxSize  = 500)
                
                sig_results   <- fgsea_results[padj < padj_thresh]
                top_up        <- head(sig_results[order(-NES)]$pathway, 10)
                top_down      <- head(sig_results[order(NES)]$pathway, 10)
                top_pathways  <- unique(c(top_up, top_down))
                
                if (!is_empty(top_pathways)) {
                        p <- as.ggplot(plotGseaTable(cytopus_list[top_pathways],
                                                     rank, 
                                                     fgsea_results,
                                                     gseaParam = 0.5)) +
                                labs(title = paste0("GSEA of DEGs according to ",
                                                    cfg$title_desc,
                                                    " for ",
                                                    cell_type,
                                                    " cells using cytopus gene set."))
                        
                        # Define output directory and filename prefix
                        out_dir     <- paste0(out_dir, cfg$out_folder, "/")
                        file_prefix <- paste0(out_dir, gsub("/", "_", cell_type))
                        
                        ggsave(p, filename = paste0(file_prefix, ".pdf"), device = "pdf", width = 12, height = 8, create.dir = TRUE)
                        write_csv(sig_results, file = paste0(file_prefix, ".csv"))
                        
                        # clusterProfiler-style dotplot
                        plot_data <- sig_results[pathway %in% top_pathways]
                        
                        # Convert pathway to a factor ordered by NES
                        plot_data$pathway <- factor(plot_data$pathway, levels = plot_data$pathway[order(plot_data$NES)])
                        
                        p_cp_mimic <- ggplot(plot_data, aes(x = NES, y = pathway, color = padj, size = size)) +
                                geom_point() +
                                scale_color_gradient(low = "red", high = "blue", name = "p.adjust", guide = guide_colorbar(reverse = TRUE)) +
                                scale_size_continuous(range = c(3, 8), name = "Set Size") +
                                theme_minimal() +
                                labs(title = paste0("GSEA of DEGs according to ", cfg$title_desc, "\nfor ", cell_type, " cells using cytopus gene set."),
                                     x = "Normalized Enrichment Score (NES)",
                                     y = NULL) +
                                theme(axis.text.y = element_text(size = 10))
                        
                        ggsave(p_cp_mimic, filename = paste0(file_prefix, "_dotplot.pdf"), device = "pdf", width = 12, height = 8, create.dir = TRUE)
                }
        }
}
write("Done.", file = snakemake@output[[1]])
message("Done")

# Close Logging ----------------------------------------------------------------
sink(type="message")
sink()