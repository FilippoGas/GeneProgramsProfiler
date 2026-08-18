# ==============================================================================
# Script: spectra_WMW.R
# Author: Filippo Gastaldello
# Date: 06/08/2026
# Description: 
#   Analyze spectra output and looks for deregulated gene programs. 
#   Gene program activations are tested for differences between the condition,
#   per each celltype and per each patient using Wilcoxon-Mann-Whitney U-test
#
# Snakemake Expected Inputs:
#   - snakemake@input[["cell_scores"]] : Path to the cell score matrix
#                                        (spectra's output)
# Snakemake Expected Outputs:
#   - snakemake@output[["gp_activation_WMW"]] : Results from the WMW test
#                                       
# Snakemake Expected Params:
#   - snakemake@params[["case"]] : Name of the case condition in the metadata
#   - snakemake@params[["control]] : Name of the control condition in the 
#                                    metadata
# ==============================================================================

# Setup Logging ----------------------------------------------------------------
# Redirect all output and messages to the Snakemake log file
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# Load Libraries ---------------------------------------------------------------
suppressPackageStartupMessages({
        library(tidyverse)
        library(parallel)
})

message("Starting R script \"spectra_WMW.R\"...")

# 1. Load Data ####
# ------------------------------------------------------------------------------
message("Loading cell score matrix from: ", snakemake@input[["cell_scores"]])
cell_scores <- read_csv(snakemake@input[["cell_scores"]])
message("Done")
message("Loading case condition name: ", snakemake@params[["case"]])
case <- snakemake@params[["case"]]
message("Done")
message("Loading control condition name: ", snakemake@params[["control"]])
control <- snakemake@params[["control"]]
message("Done")

# 2. Test gene program activation differences between case and coltrol ####
# ------------------------------------------------------------------------------
message("Starting Wilcoxon-Mann-Whitney U-test to look for deregulated
        gene programs ...")
res_WMW <- mclapply(unique(cell_scores$celltype), function(cell_type){
        
        # Subset cell scores to keep only current cell type
        celltype_cell_score <- cell_scores %>%
                dplyr::filter(celltype == cell_type)
        
        # Store total cells for percentages
        n_total_control <- sum(celltype_cell_score$diagnosis == control)
        n_total_ipf     <- sum(celltype_cell_score$diagnosis == case)
        
        # Use a list to store loop results
        tests_list <- list()
        
        programs <- colnames(
                celltype_cell_score)[6:length(colnames(celltype_cell_score))]
        
        for (program in programs) {
                
                # Extract active cells
                vec_control <- celltype_cell_score %>%
                        dplyr::filter(diagnosis == "Control" & 
                                      .[[program]] > 0.001) %>%
                        pull(all_of(program))
                vec_ipf     <- celltype_cell_score %>%
                        dplyr::filter(diagnosis == "IPF" &
                                      .[[program]] > 0.001) %>%
                        pull(all_of(program))
                
                n_active_control <- length(vec_control)
                n_active_case    <- length(vec_ipf)
                
                perc_active_control <- round((n_active_control / 
                                              n_total_control) * 100, 2)
                perc_active_case    <- round((n_active_case /
                                              n_total_ipf) * 100, 2)
                
                # Try to compute the test
                test_res <- try(wilcox.test(x = vec_control,
                                            y = vec_ipf,
                                            paired = FALSE),
                                silent = TRUE)
                
                if (inherits(test_res, "try-error")) {
                        stat <- NA
                        pval <- NA
                }else{
                        stat <- test_res$statistic
                        pval <- test_res$p.value
                }
                
                # Store results, including the active sample sizes
                tests_list[[program]] <- data.frame(
                        cell_type = cell_type, 
                        program = program, 
                        active_cells_control_perc = perc_active_control,
                        active_cells_case_perc = perc_active_case, 
                        n_active_control = n_active_control, 
                        n_active_case = n_active_case,       
                        statistic = stat, 
                        p = pval
                )
        }
        return(bind_rows(tests_list))
},
mc.cores=snakemake@threads,
mc.preschedule = TRUE,
mc.cleanup=TRUE)

# Bind all cell types together
pg_activation_WMW <- bind_rows(res_WMW)

# Compute FDR and Effect Size
pg_activation_WMW <- pg_activation_WMW %>%
        drop_na(p) %>% 
        group_by(cell_type) %>% 
        mutate(FDR = p.adjust(p, method = "fdr")) %>% 
        ungroup() %>% 
        relocate(FDR, .after = p) %>%
        # Apply the Rank-Biserial Formula using the sample sizes
        mutate(effect_size = abs(1 - (2 * statistic) /
                                (n_active_control * n_active_case)))

write_csv(pg_activation_WMW, snakemake@output[["gp_activation_WMW"]])
message("Done")

# Close Logging ----------------------------------------------------------------
sink(type="message")
sink()