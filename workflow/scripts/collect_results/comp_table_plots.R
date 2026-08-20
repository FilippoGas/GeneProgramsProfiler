# ==============================================================================
# Script: comp_table_plots.R
# Author: Filippo Gastaldello
# Date: 20/08/2026
# Description: 
#   Make concordance plots showing detected gene program deregulation
#   across methods.
#
# Snakemake Expected Inputs:
#   - snakemake@input[["table"]]: Path to comparison table
#
# Snakemake Expected Outputs:
#   - snakemake@output[["horizontal_upset"]] : Path to save horizontal upset 
#                                              plot on methods concordance
#   - snakemake@output[["vertical_upset"]] : Path to save vertical upset 
#                                            plot on methods concordance
#   - snakemake@output[["intersections"]] : Path to save the intersections csv, 
#                                           containing the full list of programs
#                                           in each intersection
#
# Snakemake Expected Params:
#   - snakemake@params[["FDR_thresh"]]: FDR threshold for spectra and cNMF results
#   - snakemake@params[["log2FC_thresh"]] : log2FC threshold for spectra and cNMF results
#   - snakemake@params[["effect_size_thresh"]] : effect size threshold for spectra and cNMF results
#                                       
# ==============================================================================

# Setup Logging ----------------------------------------------------------------
# Redirect all output and messages to the Snakemake log file
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# Prevent R from generating the default Rplots.pdf file
pdf(NULL)

# Load Libraries ---------------------------------------------------------------
suppressPackageStartupMessages({
        library(tidyverse)
        library(UpSetR)
        library(ComplexHeatmap)
        library(openxlsx)
})

message("Starting R script \"comp_table_plots.R\"...")

# 1. Load Data ####
# ------------------------------------------------------------------------------
message("Loading comparison table from: ", snakemake@input[["table"]])
comp_table <- read_csv(snakemake@input[["table"]])
message("Done")
message("Loading FDR threshold for Wilcoxon-Mann-whitney and Linea Mixed Models results: ", snakemake@params[["FDR_thresh"]])
FDR_thresh <- snakemake@params[["FDR_thresh"]]
message("Done")
message("Loading log2FC threshold for Linea Mixed Models results: ", snakemake@params[["log2FC_thresh"]])
log2FC_thresh <- snakemake@params[["log2FC_thresh"]]
message("Done")
message("Loading effect size threshold for Wilcoxon-Mann-whitney results: ", snakemake@params[["effect_size_thresh"]])
effect_size_thresh <- snakemake@params[["effect_size_thresh"]]
message("Done")

# 3. UpSet plot of methods concordance ####
# ------------------------------------------------------------------------------
message("Creating concordance list ...")
upset_data <- list(
        spectra_WMW = comp_table %>%
                filter(FDR_wmw_spectra < FDR_thresh,
                       abs(effect_size_wmw_spectra) > effect_size_thresh) %>% 
                mutate(label = paste0(cell_type, "_", program_short)) %>% 
                pull(label) %>% 
                unique(),
        
        spectra_LMM = comp_table %>% 
                filter(FDR_lmm_spectra<FDR_thresh,
                       abs(log2FC_lmm_spectra)>log2FC_thresh) %>% 
                mutate(label = paste0(cell_type, "_", program_short)) %>% 
                pull(label) %>% 
                unique(),
        
        cNMF_WMW = comp_table %>%
                filter(FDR_wmw_cNMF < FDR_thresh,
                       abs(effect_size_wmw_cNMF) > effect_size_thresh) %>% 
                mutate(label = paste0(cell_type, "_", program_short)) %>% 
                pull(label) %>% 
                unique(),
        
        cNMF_LMM = comp_table %>% 
                filter(FDR_lmm_cNMF<FDR_thresh,
                       abs(log2FC_lmm_cNMF)>log2FC_thresh) %>% 
                mutate(label = paste0(cell_type, "_", program_short)) %>% 
                pull(label) %>% 
                unique(),
        
        GSEA_sc = comp_table %>% 
                filter(!is.na(padj_GSEA_sc)) %>% 
                mutate(label = paste0(cell_type, "_", program_short)) %>% 
                pull(label) %>% 
                unique(),
        
        GSEA_pb = comp_table %>%
                filter(!is.na(padj_GSEA_pb)) %>%
                mutate(label = paste0(cell_type, "_", program_short)) %>% 
                pull(label) %>% 
                unique(),
        
        ORA_sc = comp_table %>% 
                filter(!is.na(padj_ORA_UP_sc) | !is.na(padj_ORA_DOWN_sc)) %>% 
                mutate(label = paste0(cell_type, "_", program_short)) %>% 
                pull(label) %>% 
                unique(),
        
        ORA_pb = comp_table %>%
                filter(!is.na(padj_ORA_UP_pb) | !is.na(padj_ORA_DOWN_pb)) %>% 
                mutate(label = paste0(cell_type, "_", program_short)) %>% 
                pull(label) %>% 
                unique(),
        
        ORA_both = comp_table %>% 
                filter(!is.na(padj_ORA_UP_both) | !is.na(padj_ORA_DOWN_both)) %>% 
                mutate(label = paste0(cell_type, "_", program_short)) %>% 
                pull(label) %>% 
                unique()
)
message("Done")
message("Plotting horizontal UpSet plot ...")
# Filter out completely empty sets to prevent UpSetR from crashing
upset_data_filtered <- upset_data[lengths(upset_data) > 0]

# Only attempt to plot if at least 2 methods found significant programs
if(length(upset_data_filtered) >= 2) {
        
        df_upset <- fromList(upset_data_filtered)
        
        pdf(file = snakemake@output[["horizontal_upset"]], width = 9, height = 5, onefile = FALSE)
        
        # 'sets' parameter forces UpSetR to include all active columns, 
        # preventing it from silently dropping smaller sets like LMM from the visualization.
        print(upset(
                df_upset,
                order.by = "freq",
                sets = names(upset_data_filtered), 
                nsets = length(upset_data_filtered)
        ))
        
        dev.off()
        
} else {
        message("Warning: Fewer than 2 methods yielded significant programs. Outputting empty placeholder plot.")
        pdf(file = snakemake@output[["horizontal_upset"]], width = 9, height = 5, onefile = FALSE)
        plot.new()
        text(0.5, 0.5, "Not enough overlapping methods to compute UpSet plot")
        dev.off()
}
message("Done")

message("Plotting vertical UpSet plot ...")
# Define the exact order you want for the columns (sets)
ordered_sets <- c(
        "ORA_sc", "ORA_pb", "ORA_both",         # ORA methods
        "GSEA_sc", "GSEA_pb",                   # GSEA methods
        "spectra_WMW", "cNMF_WMW",              # WMW methods
        "spectra_LMM", "cNMF_LMM"               # LMM methods
)

# Create combination matrix keeping all intersections
m <- make_comb_mat(upset_data)
m_transposed <- t(m)

all_codes <- comb_name(m_transposed)
all_sizes <- comb_size(m_transposed)
all_degrees <- comb_degree(m_transposed) 
plot_order <- order(all_degrees, all_sizes, decreasing = TRUE)

# Map sequential indices so Row 1 in the plot = Index 1 in the table
idx_labels <- integer(length(plot_order))
idx_labels[plot_order] <- 1:length(plot_order)
plot_index_strings <- paste0("[ ", idx_labels, " ]")

# Build the text labels and the export table (Intersection_Code removed)
intersection_table <- data.frame(
        ID = integer(),
        Degree = numeric(),
        Size = numeric(),
        Items = character()
)

all_labels <- sapply(seq_along(all_codes), function(i) {
        items <- extract_comb(m, all_codes[i])
        
        # Save full data to the export table
        intersection_table <<- bind_rows(
                intersection_table,
                data.frame(
                        ID = idx_labels[i],
                        Degree = all_degrees[i],
                        Size = all_sizes[i],
                        Items = paste(items, collapse = "\n")
                )
        )
        
        # Create truncated labels for the plot
        if (length(items) == 0) {
                return("")
        } else if (length(items) <= 2) {
                return(paste(items, collapse = "\n"))
        } else {
                truncated <- paste(items[1:2], collapse = "\n")
                return(paste0(truncated, "\n(+ ", length(items) - 2, " more)"))
        }
})

# Sort the table so it matches the plot order top-to-bottom
intersection_table <- intersection_table %>% arrange(ID)

# Create the Annotations
custom_left_annotation <- rowAnnotation(
        ID = anno_text(
                plot_index_strings,
                just = "left",                 # Strictly align text to the left
                location = unit(5, "mm"),      # Guarantee exactly 5mm of blank space on the left
                gp = gpar(fontsize = 10, fontface = "bold")
        ),
        annotation_label = "ID",               # Explicitly define the title
        show_annotation_name = TRUE,           # Force the title to show
        annotation_name_side = "top",          # Place it at the top
        annotation_name_rot = 0,               # Keep it horizontal
        annotation_name_gp = gpar(fontface = "bold", fontsize = 12),
        width = unit(2.0, "cm")                # Widen the bounding box to ensure no clipping
)

custom_right_annotation <- rowAnnotation(
        Size = anno_text(
                as.character(all_sizes), 
                just = "center",
                gp = gpar(fontsize = 10, fontface = "bold")
        ),
        Programs = anno_text(
                all_labels,
                just = "left",
                location = unit(2, "mm"),
                gp = gpar(fontsize = 7) 
        ),
        show_annotation_name = TRUE,
        annotation_name_side = "top", 
        annotation_name_rot = 0,
        annotation_name_gp = gpar(fontface = "bold", fontsize = 12)
)

# Calculate dynamic dimensions and draw the plot
n_intersections <- length(plot_order)
dynamic_height <- max(8, 5 + (n_intersections * 0.4)) 
n_methods <- length(ordered_sets)
dynamic_width <- max(10, 8 + (n_methods * 0.5)) 

message("Dynamic plot dimensions calculated: ", round(dynamic_width, 1), "W x ", round(dynamic_height, 1), "H inches")

pdf(file = snakemake@output[["vertical_upset"]], width = dynamic_width, height = dynamic_height, onefile = FALSE)

# Adding padding = unit(c(bottom, left, top, right)) prevents top/left labels from being clipped by the PDF edge
draw(UpSet(
        m_transposed, 
        set_order = ordered_sets,               
        comb_order = plot_order,
        left_annotation = custom_left_annotation,
        right_annotation = custom_right_annotation,
        column_names_side = "top",
        column_names_rot  = 45
), padding = unit(c(5, 5, 15, 5), "mm")) 

dev.off()
message("Done")

# Export as a formatted Excel table with Text Wrapping
message("Saving intersections Excel table to: ", snakemake@output[["intersections"]])

wb <- createWorkbook()
addWorksheet(wb, "Intersections")
writeData(wb, "Intersections", intersection_table)

# Apply text wrap styling (Updated to cols = 4 since Code column is gone)
wrap_style <- createStyle(wrapText = TRUE, valign = "top")
addStyle(wb, "Intersections", style = wrap_style, 
         rows = 1:(nrow(intersection_table) + 1), 
         cols = 4, 
         gridExpand = TRUE)

# Widen the Items column to make it readable
setColWidths(wb, "Intersections", cols = 4, widths = 50)

saveWorkbook(wb, snakemake@output[["intersections"]], overwrite = TRUE)
message("Done")

message("Saving intersections table to: ", snakemake@output[["intersections"]])
write_csv(intersection_table, snakemake@output[["intersections"]])
message("Done")

# Close Logging ----------------------------------------------------------------
sink(type="message")
sink()