# ==============================================================================
# Script: DEA_plots.R
# Author: Filippo Gastaldello
# Date: 04/08/2026
# Description: 
#   Plots some statistics and summaries on the results of the DE analysis
#
# Snakemake Expected Inputs:
#   - snakemake@input[["DEGs_sc"]]   : Path to results of single cell DE 
#                                      analysis
#   - snakemake@input[["DEGs_pb"]]   : Path to results of pseudo-bulk DE
#                                      analysis
#   - snakemake@input[["DEGs_both"]] : Path to merged results (sc and pb)
#
# Snakemake Expected Outputs:
#   - snakemake@output[["overlap]] : Out path for single cell/pseudobulk 
#                                    overlap plot 
#   - snakemake@output[["FDR_correlation.pdf"]] : Out path for correlation
#                                                 between FDR in single cell vs
#                                                 pseudo bulk
#   - snakemake@output[["log2FC_correlation.pdf"]] : Out path for correlation
#                                                    between log2FC in single
#                                                    cell vs pseudo bulk
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
        library(eulerr)
        library(ggplotify)
        library(ggrastr)
        library(patchwork)
})

message("Starting R script \"DEA_plots.R\"...")

# 1. Load Data ####
# ------------------------------------------------------------------------------
message("Loading results from single cell differential analysis from ",
        snakemake@input[["DEGs_sc"]])
DEGs_sc <- read_csv(snakemake@input[["DEGs_sc"]])
message("Done")
message("Loading results from pseudo-bulk differential analysis from ",
        snakemake@input[["DEGs_pb"]])
DEGs_pb <- read_csv(snakemake@input[["DEGs_pb"]])
message("Done")
message("Loading merged results from ", snakemake@input[["DEGs_both"]])
DEGs_both <- read_csv(snakemake@input[["DEGs_both"]])
message("Done")

# 2. Check overlap between single cell and pseudo-bulk analysis ####
# ------------------------------------------------------------------------------
message("Plotting overlap between single cell and pseudo bulk analysis")
plot_list <- lapply(unique(c(DEGs_sc$celltype, DEGs_pb$celltype)), 
                    function(cell_type) {
        sc_genes <- DEGs_sc %>% 
                        dplyr::filter(is_sig == TRUE,
                                      celltype == cell_type) %>%
                        pull(gene)
        pb_genes <- DEGs_pb %>%
                        dplyr::filter(is_sig == TRUE,
                                      celltype == cell_type) %>%
                        pull(gene)
        fit <- euler(list("sc" = sc_genes,
                          "Pseudo-bulk" = pb_genes)
                     )
        # Plot ven
        p1 <- as.ggplot(plot(fit,
                             quantities = TRUE,
                             fills = c("skyblue", "pink"),
                             alpha = 0.5)
                        )
        # Combine with correlation of FC in sc vs pb
        p2 <- DEGs_both %>%
                dplyr::filter(both_sig,
                              celltype == cell_type) %>% 
                dplyr::select(avg_log2FC_sc,
                              avg_log2FC_pb) %>%
                drop_na() %>%
                ggplot(aes(x = avg_log2FC_sc, y = avg_log2FC_pb)) +
                geom_point_rast(
                        size = 0.1
                        ) +
                geom_smooth(
                        method=lm,
                        color="red",
                        fill="#69b3a2",
                        se=TRUE,
                        linewidth = 0.25
                        ) +
                labs(
                        title = cell_type
                        ) +
                theme_minimal()
        p <- wrap_plots(c(p1, p2)) + plot_annotation(title = cell_type)
})

combined_plot <- wrap_plots(plot_list) + 
        plot_annotation(title    = "Overlap between DEGs_sc found in scDE 
                                    and pseudo-bulkDE analyses",
                        subtitle = "Scatterplot showing logFC correlation
                                    between genes significantly altered
                                    (FDR<0.05 and log2FC>1) in both analyses
                                    (venn intersection)"
                        )
ggsave(filename = snakemake@output[["overlap"]],
       plot = combined_plot,
       width = 25,
       height = 11,
       create.dir = TRUE)
message("Done")

# 3. Plot correlation between FDR from single cell and pseudo bulk
# ------------------------------------------------------------------------------
message("Plotting correlation between FDR in single cell and pseudo bulk")
# Significant in one of the 2
p1 <- DEGs_both %>%
        dplyr::filter(is_sig_sc | is_sig_pb) %>% 
        dplyr::select(FDR_sc, FDR_pb) %>%
        drop_na() %>%
        dplyr::filter(!FDR_sc == 0) %>% 
        ggplot(aes(x = -log(FDR_sc), y = -log(FDR_pb))) +
                geom_point_rast(
                        size = 0.01,
                        alpha = 0.07,
                        raster.dpi = 400
                        ) +
                geom_smooth(
                        method=lm ,
                        color="red",
                        fill="#69b3a2",
                        se=TRUE,
                        size = 0.25
                        ) +
                labs(
                        title = "Significant in at least one of the two
                                 analyses"
                        ) +
                theme_minimal()

# Significant in both
p2 <- DEGs_both %>%
        dplyr::filter(both_sig) %>% 
        dplyr::select(FDR_sc, FDR_pb) %>%
        drop_na() %>%
        dplyr::filter(!FDR_sc == 0) %>% 
        ggplot(aes(x = -log(FDR_sc), y = -log(FDR_pb))) +
                geom_point_rast(
                        size = 0.01,
                        alpha = 0.07,
                        raster.dpi = 400
                        ) +
                geom_smooth(
                        method=lm ,
                        color="red",
                        fill="#69b3a2",
                        se=TRUE,
                        size = 0.25
                        ) +
                labs(
                        title = "Significant in both analyses"
                        ) +
        theme_minimal()
p <- wrap_plots(p1, p2) +
        plot_annotation(
                title = "Correlation between single cell and pseudobulk
                         FDRs of differential analysis",
                subtitle = "(Removed genes where single cell FDR was zero)")
ggsave(p,
       filename = snakemake@output[["FDR_correlation.pdf"]],
       device = "pdf",
       width = 10,
       height = 5,
       create.dir = TRUE)
message("Done")

# 4. Plot correlation between log2FC from single cell and pseudo bulk
# ------------------------------------------------------------------------------
message("Plotting correlation between log2FC in single cell and pseudo bulk")
# Significant in one of the 2
p1 <- DEGs_both %>%
        dplyr::filter(is_sig_sc | is_sig_pb) %>% 
        dplyr::select(avg_log2FC_sc, avg_log2FC_pb) %>%
        drop_na() %>% 
        ggplot(aes(x = avg_log2FC_sc, y = avg_log2FC_pb)) +
                geom_point_rast(
                        size = 0.01,
                        alpha = 0.07,
                        raster.dpi = 400
                        ) +
                geom_smooth(
                        method=lm,
                        color="red",
                        fill="#69b3a2",
                        se=TRUE,
                        size = 0.25
                        ) +
                labs(
                        title = "Significant in at least one of the two 
                                 analyses"
                        ) +
                theme_minimal()
# Significant in both
p2 <- DEGs_both %>%
        dplyr::filter(both_sig) %>% 
        dplyr::select(avg_log2FC_sc, avg_log2FC_pb) %>%
        drop_na() %>%
        ggplot(aes(x = avg_log2FC_sc, y = avg_log2FC_pb)) +
                geom_point_rast(
                        size = 0.01,
                        alpha = 0.07,
                        raster.dpi = 400
                        ) +
                geom_smooth(
                        method=lm,
                        color="red",
                        fill="#69b3a2",
                        se=TRUE,
                        size = 0.25
                        ) +
                labs(
                        title = "Significant in both analyses"
                        ) +
                theme_minimal()
p <- wrap_plots(p1, p2) +
        plot_annotation(
                title = "Correlation between single cell and pseudobulk
                         Log2FC of differential analysis"
                )
ggsave(p,
       filename = snakemake@output[["log2FC_correlation.pdf"]],
       device = "pdf",
       width = 10,
       height = 5,
       create.dir = TRUE)
message("Done")

# Close Logging ----------------------------------------------------------------
sink(type="message")
sink()