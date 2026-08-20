# ==============================================================================
# Script: WMW.R
# Author: Filippo Gastaldello
# Date: 06/08/2026
# Description: 
#   Analyze gene program activation data and look for deregulated gene programs. 
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
#   - snakemake@params[["active_cell_thresh"]] : Activation threshold to 
#                                                consider a program active in
#                                                cell
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
        library(effectsize)
})

message("Starting R script \"WMW.R\"...")

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
message("Loading program activation threshold: ",
        snakemake@params[["active_cell_thresh"]])
active_cell_thresh <- snakemake@params[["active_cell_thresh"]]
message("Done")

# 2. Test gene program activation differences between case and coltrol ####
# ------------------------------------------------------------------------------
message("Starting Wilcoxon-Mann-Whitney U-test to look for deregulated
        gene programs ...")
res_WMW_sample_agg <- lapply(unique(cell_scores$celltype), function(cell_type){
        
        # Subset cell scores to keep only current cell type
        celltype_cell_score <- cell_scores %>% dplyr::filter(celltype == cell_type)
        
        # Store total cells for percentages
        n_total_control <- sum(celltype_cell_score$diagnosis == control)
        n_total_ipf     <- sum(celltype_cell_score$diagnosis == case)
        
        tests_list <- list()
        
        programs <- colnames(celltype_cell_score)[6:length(colnames(celltype_cell_score))]
        
        for (program in programs) {
                
                # Extract active cells
                vec_control <- celltype_cell_score %>% 
                        dplyr::filter(diagnosis == control & .[[program]] > active_cell_thresh) %>%
                        dplyr::select(all_of(program), sample_name)
                vec_ipf     <- celltype_cell_score %>%
                        dplyr::filter(diagnosis == case & .[[program]] > active_cell_thresh) %>%
                        dplyr::select(all_of(program), sample_name)
                
                n_active_control <- nrow(vec_control)
                n_active_case    <- nrow(vec_ipf)
                
                perc_active_control <- round((n_active_control / n_total_control) * 100, 2)
                perc_active_case    <- round((n_active_case / n_total_ipf) * 100, 2)
                
                # Aggregate cell scores per sample
                vec_control_agg <- vec_control %>% summarise(.by = sample_name, score = mean(.data[[program]]))
                vec_ipf_agg     <- vec_ipf     %>% summarise(.by = sample_name, score = mean(.data[[program]]))
                
                # Extract the numeric vectors for the test
                scores_control <- vec_control_agg %>% pull(score)
                scores_ipf     <- vec_ipf_agg %>% pull(score)
                
                # Try to compute the test
                test_res <- try(wilcox.test(x          = scores_control,
                                            y          = scores_ipf,
                                            conf.int   = TRUE,
                                            paired     = FALSE),
                                silent = TRUE)
                
                if (inherits(test_res, "try-error") || length(scores_control) == 0 || length(scores_ipf) == 0) {
                        stat       <- NA
                        pval       <- NA
                        es         <- NA
                        es_ci_low  <- NA
                        es_ci_high <- NA
                } else {
                        stat <- test_res$statistic
                        pval <- test_res$p.value
                        
                        # Compute effect size and 95% CI on the sample-aggregated data
                        es_res <- effectsize::rank_biserial(x = scores_ipf, 
                                                            y = scores_control,  # The second group is used as reference by default 
                                                            ci = 0.95)
                        
                        es         <- es_res$r_rank_biserial
                        es_ci_low  <- es_res$CI_low
                        es_ci_high <- es_res$CI_high
                }
                
                # Store results
                tests_list[[program]] <- data.frame(
                        cell_type                 = cell_type, 
                        program                   = program, 
                        active_cells_control_perc = perc_active_control,
                        active_cells_case_perc    = perc_active_case, 
                        n_active_control          = n_active_control, 
                        n_active_case             = n_active_case,
                        n_samples_control         = length(scores_control),
                        n_samples_case            = length(scores_ipf),
                        avg_activation_control    = mean(scores_control),
                        avg_activation_case       = mean(scores_ipf),
                        statistic                 = stat, 
                        p                         = pval,
                        effect_size               = es,
                        effect_size_CI_lower      = es_ci_low,
                        effect_size_CI_upper      = es_ci_high
                )
        }
        return(bind_rows(tests_list))
})

# Bind all cell types together
gp_activation_WMW <- bind_rows(res_WMW_sample_agg)

# Compute FDR and Effect Size
gp_activation_WMW <- gp_activation_WMW %>%
        drop_na(p) %>% 
        group_by(cell_type) %>% 
        mutate(FDR = p.adjust(p, method = "fdr")) %>% 
        ungroup() %>% 
        relocate(FDR, .after = p)

write_csv(gp_activation_WMW, snakemake@output[["gp_activation_WMW"]])
message("Done")

# Close Logging ----------------------------------------------------------------
sink(type="message")
sink()