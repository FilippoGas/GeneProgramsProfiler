# ==============================================================================
# Script: LMM_plots.R
# Author: Filippo Gastaldello
# Date: 18/08/2026
# Description: 
#   Plot results of differential gene program activation analysis 
#
#   - snakemake@input[["gp_activation_LMM"]] : Path to differential activation
#                                              analysis results
# Snakemake Expected Outputs:
#   - snakemake@input[["volcano"]] : Path to save volcano plot of deregulated
#                                    programs
# Snakemake Expected Params:
#   - snakemake@params[["log2FC_thresh"]] : log2FC threshold for LMM's effect
#                                           size
#   - snakemake@params[["FDR_thresh]] : False Discovery Rate threshold for
#                                       Linear Mixed Models
# ==============================================================================

# Setup Logging ----------------------------------------------------------------
# Redirect all output and messages to the Snakemake log file
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# Load Libraries ---------------------------------------------------------------
suppressPackageStartupMessages({
        library(tidyverse)
        library(ggrepel)
        library(ggpubr)
})

message("Starting R script \"LMM_plots.R\"...")

# Prevent R from generating the default Rplots.pdf file
pdf(NULL)

# 1. Load Data ####
# ------------------------------------------------------------------------------
message("Loading differential activation dataset from: ",
        snakemake@input[["gp_activation_LMM"]])
GP_activation_LMM <- read_csv(snakemake@input[["gp_activation_LMM"]])
message("Done")
message("Loading log2FC threshold: ",
        snakemake@params[["log2FC_thresh"]])
log2FC_thresh <- snakemake@params[["log2FC_thresh"]]
message("Done")
message("Loading FDR threshold: ", snakemake@params[["FDR_thresh"]])
FDR_thresh <- snakemake@params[["FDR_thresh"]]
message("Done")

# 2. Plot volcano plot ####
# ------------------------------------------------------------------------------
plot_data <- GP_activation_LMM %>%
        dplyr::filter(!log2FC=="NaN",
                      !str_detect(term, "phase")) %>% 
        mutate(tag = paste0(cell_type, "-", str_split_i(program, "-X-",3)),
               diff_activated = ifelse(FDR<FDR_thresh & log2FC>log2FC_thresh, "UP", "NO"),
               diff_activated = ifelse(FDR<FDR_thresh & log2FC< -log2FC_thresh, "DOWN", diff_activated))
if(!nrow(GP_activation_LMM)==0){
        p <- plot_data %>% ggplot(aes(x=log2FC, y=-log10(FDR), color = diff_activated)) +
                geom_point(size = 1) +
                scale_color_manual(values = c("UP"="red", "NO"="grey", "DOWN"="purple"),
                                   labels = c("Downactivated", "Not significant", "Overactivated")) +
                geom_vline(xintercept = c(-log2FC_thresh, log2FC_thresh), col = 'grey', linetype = 'dashed') +
                geom_hline(yintercept = -log10(FDR_thresh), col = 'grey', linetype = 'dashed') +
                labs(
                        title = "Linear Mixed Model - Sample ID as random effect, correcting for cell cycle phase",
                        subtitle = "comparison between activation in case vs control, for each program in each cell type",
                        x = "log2(activation FC)"
                ) +
                geom_label_repel(
                        data = plot_data %>% filter((FDR<FDR_thresh & log2FC>log2FC_thresh) | (FDR<FDR_thresh & log2FC< -log2FC_thresh)),
                        aes(label = tag),
                        size = 2.5,           
                        max.overlaps = Inf,   
                        box.padding = 0.5,
                        point.padding = 0.2,
                        min.segment.length = 0,
                        nudge_x = .15,
                        nudge_y = .5,) +
                theme_minimal()
        p
        ggsave(filename = snakemake@output[["volcano"]],plot = p, device = "pdf", width = 10, height = 10, create.dir = TRUE)
}else{
        message("No deregulated programs found by LMM")
        write(x="No program was found deregulated.", file = snakemake@output[["volcano"]])
}
# Close Logging ----------------------------------------------------------------
sink(type="message")
sink()