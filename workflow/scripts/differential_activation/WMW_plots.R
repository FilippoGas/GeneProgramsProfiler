# ==============================================================================
# Script: WMW_plots.R
# Author: Filippo Gastaldello
# Date: 07/08/2026
# Description: 
#   Plot results of differential gene program activation analysis 
#
# Snakemake Expected Inputs:
#   - snakemake@input[["cell_scores]] : Path to cell scores matrix
#   - snakemake@input[["gp_activation_WMW"]] : Path to differential activation
#                                              analysis results
#   - snakemake@input[["output"]] : Path to output placeholder
# Snakemake Expected Params:
#   - snakemake@params[["effect_size_thresh"]] : "Effect size threshold for 
#                                                 Wilcoxon-Mann-Whitney U-test's
#                                                 rank-biserial correlation"
#   - snakemake@params[["FDR_thresh]] : "False Discovery Rate threshold for
#                                        Wilcoxon-Mann-Whitney U-test"
# ==============================================================================

# Setup Logging ----------------------------------------------------------------
# Redirect all output and messages to the Snakemake log file
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# Load Libraries ---------------------------------------------------------------
suppressPackageStartupMessages({
        library(tidyverse)
        library(ggpubr)
        library(ggrepel)
})

message("Starting R script \"WMW_plots.R\"...")

# Prevent R from generating the default Rplots.pdf file
pdf(NULL)

# 1. Load Data ####
# ------------------------------------------------------------------------------
message("Loading cell scores matrix from: ",
        snakemake@input[["cell_scores"]])
cell_scores <- read_csv(snakemake@input[["cell_scores"]])
message("Done")
message("Loading differential activation dataset from: ",
        snakemake@input[["gp_activation_WMW"]])
gp_activation_WMW <- read_csv(snakemake@input[["gp_activation_WMW"]])
message("Done")
message("Loading effect size threshold: ",
        snakemake@params[["effect_size_thresh"]])
effect_size_thresh <- snakemake@params[["effect_size_thresh"]]
message("Done")
message("Loading FDR threshold: ", snakemake@params[["FDR_thresh"]])
FDR_thresh <- snakemake@params[["FDR_thresh"]]
message("Done")
message("Loading base output dir :",
        str_split_i(snakemake@output[[1]], ".done",1))
out_dir <- str_split_i(snakemake@output[[1]], ".done",1)
message("Done")

if(!nrow(gp_activation_WMW)==0){
        # 2. PLot combined boxplot of deregulated programs for all cell types ####
        # ------------------------------------------------------------------------------
        message("Plotting combined boxplot for all cell types ...")
        # Extract significant tests
        WMW_sample_agg_sig <- gp_activation_WMW %>%
                dplyr::filter(FDR < FDR_thresh & 
                              abs(effect_size) > effect_size_thresh)
        
        plot_data <- cell_scores %>%
                pivot_longer(6:length(colnames(cell_scores)), names_to = "program") %>%                           
                filter(value > 0.001) %>% # Only keep active cells
                group_by(sample_name, celltype, program, diagnosis) %>% 
                summarize(value = mean(value)) %>%                                                                
                ungroup() %>% 
                inner_join(WMW_sample_agg_sig %>%                                                                 
                                   dplyr::select(cell_type, program, FDR),
                           by = join_by("celltype" == "cell_type",
                                        "program" == "program"))
        
        p <- plot_data %>% ggplot(aes(x=program, y=value, fill=diagnosis))+
                geom_boxplot() +
                theme_minimal() +
                theme(
                        axis.text.x = element_text(hjust=1, angle = 45, size = 5)
                ) +
                facet_wrap(~celltype, scale="free")
        
        ggsave(filename = paste0(out_dir, "boxplot_combined.pdf"),
               plot = p,
               device = "pdf",
               width = 18,
               height = 15,
               create.dir = TRUE)
        message("Done")
        
        # 3. Plot boxplots of deregulated programs for each cell type ####
        # ------------------------------------------------------------------------------
        message("Plotting boxplots for each cell type ...")
        for (cell_type in unique(plot_data$celltype)) {
                
                p <- plot_data %>%
                        dplyr::filter(celltype==cell_type) %>%
                        ggplot(aes(x=program, y=value, fill=diagnosis)) +
                        geom_boxplot(
                                outliers = FALSE
                        ) +
                        geom_jitter(
                                position = position_jitterdodge(),
                                alpha    = 0.5,
                                size     = 0.7
                        ) +
                        labs(
                                title = paste0("IGP with significant differences in activation between case and control in ",
                                               cell_type," cells."),
                                y     = "Activation level"
                        ) +
                        stat_compare_means(
                                aes(group = diagnosis),
                                method = "wilcox.test",
                                method.args = list(paired=FALSE),
                                vjust = -1,
                                label = "p.signif"
                        ) +
                        stat_summary(
                                fun.data = function(y) {
                                        return(data.frame(
                                                y = max(y),
                                                label = paste0("n= ",length(y))
                                        ))
                                },
                                geom = "text",
                                position = position_dodge(width = 0.75),
                                vjust = -0.5,
                                size = 3
                        ) +
                        theme_minimal() +
                        theme(
                                axis.text.x = element_text(hjust=1, angle=45, size=8)
                        )
                
                
                ggsave(filename = paste0(out_dir, gsub("/", "_", cell_type),".pdf"),
                       plot = p,
                       device = "pdf",
                       width = 12,
                       height = 8,
                       create.dir = TRUE)
        }
        message("Done")
        
        # 4. Plot volcano plot of FDR vs effet size ####
        # ------------------------------------------------------------------------------
        message("Plotting volcano plot of FDR vs effect size ...")
        plot_data <- gp_activation_WMW %>% 
                mutate(label = paste0("(", cell_type, ") ",
                                      str_split_i(program, "-X-", 3)),
                       sig = case_when(
                               FDR <= 0.05 & effect_size >= 0.5                      ~ "up_high",
                               FDR <= 0.05 & effect_size >=0.3 & effect_size < 0.5   ~ "up_mid",
                               FDR <= 0.05 & effect_size <= -0.5                     ~ "down_high",
                               FDR <= 0.05 & effect_size <= -0.3 & effect_size> -0.5 ~ "down_mid",
                               TRUE                                    ~ "no_sig"
                       )) %>%
                dplyr::select(label,
                              FDR,
                              effect_size,
                              sig)
        
        p <- plot_data %>% ggplot(aes(x = effect_size, y = -log10(FDR), color = sig)) +
                geom_point() +
                scale_color_manual(
                        values = c("up_high" = "#f00514",
                                   "up_mid" = "#d19f9f",
                                   "down_high" = "#6e02e0",
                                   "down_mid" = "#b6a9c9",
                                   "no_sig" = "#d4d4d4"),
                        breaks = c("down_high", "down_mid", "up_high", "up_mid"),
                        labels = c("Over activated - High confidence",
                                   "Over activated - Moderate confidence",
                                   "Repressed - High confidence",
                                   "Repressed - Moderate confidence")
                ) +
                geom_label_repel(
                        data = plot_data %>% 
                                filter(sig == "up_high" | sig == "down_high"),
                        aes(label = label),
                        size = 2.5,           
                        max.overlaps = Inf,   
                        box.padding = 0.5,
                        point.padding = 0.2,
                        min.segment.length = 0,
                        nudge_x = .15,
                        nudge_y = .5,
                ) +
                labs(
                        title    = "Deregulated programs activations",
                        subtitle = "Wilcoxon unpaired test on sample aggregated scores per cell type.",
                        color    = "Deregulation effect"
                ) +
                theme_minimal()
        ggsave(plot = p,
               filename = paste0(out_dir, "volcano_FDR_effect_size.pdf"),
               create.dir = TRUE,
               height = 9,
               width = 13,
               device = "pdf")
        message("Done")
}else{
        message("No deregulated program found by WMW")
}
write("Done.", file = snakemake@output[["done"]])
# Close Logging ----------------------------------------------------------------
sink(type="message")
sink()