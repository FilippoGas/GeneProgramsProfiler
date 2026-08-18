# ==============================================================================
# Script: spectra_LMM.R
# Author: Filippo Gastaldello
# Date: 18/08/2026
# Description: 
#   Analyze spectra output and looks for deregulated gene programs. 
#   Gene program activations are tested for differences between the condition,
#   per each celltype and per each patient using Linear Mixed models to correct
#   for cell cycle phase
#
# Snakemake Expected Inputs:
#   - snakemake@input[["cell_scores"]] : Path to the cell score matrix
#                                        (spectra's output)
# Snakemake Expected Outputs:
#   - snakemake@output[["gp_activation_LMM"]] : Results from the LMM test
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
        library(lme4)
        library(lmerTest)
})

message("Starting R script \"spectra_LMM.R\"...")

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

# 2. Test differential activation using LMMs ####
# ------------------------------------------------------------------------------
message("Testing differential gene program activation with Linear Mixed Models...")
res_lmm <- mclapply(unique(cell_scores$celltype),
                  function(cell_type){
                          # Subset cell scores to keep only current cell type
                          celltype_cell_score <- cell_scores %>% 
                                 dplyr::filter(celltype == cell_type)
                          # Prepare dataset to store LMM's test results
                          tests <- data_frame("cell_type"=character(),
                                              "program"=character(),
                                              "term"=character(),
                                              "active_cells_control_perc"=numeric(),
                                              "active_cells_case_perc"=numeric(),
                                              "avg_control_activation"=numeric(),
                                              "avg_case_activation"= numeric(),
                                              "estimate"=numeric(),
                                              "std.err"=numeric(),
                                              "statistic"=numeric(),
                                              "p"=numeric())
                          # For each program, test for difference in activation between case and control, only in active cells (activation > 0.001) as per Authors suggestion
                          for (program in colnames(celltype_cell_score)[6:length(colnames(celltype_cell_score))]) {
                                  # Compute number of active cells (activation > 0.001) in this cell type
                                  perc_active_cells_control <- round((celltype_cell_score %>%
                                                                              dplyr::filter(diagnosis == control & .[[program]]>0.1) %>%
                                                                              rownames() %>%
                                                                              length())/length(rownames(celltype_cell_score))*100,
                                                                     2)
                                  perc_active_cells_case <-    round((celltype_cell_score %>%
                                                                              dplyr::filter(Diagnosis == case & .[[program]]>0.1) %>%
                                                                              rownames() %>%
                                                                              length())/length(rownames(celltype_cell_score))*100,
                                                                     2)
                                  # Compute average activation in case vs control
                                  avg_control_activation <- celltype_cell_score %>%
                                          dplyr::filter(Diagnosis == control) %>% 
                                          pull(any_of(program)) %>% 
                                          mean()
                                  avg_case_activation <- celltype_cell_score %>%
                                          dplyr::filter(Diagnosis == case) %>% 
                                          pull(any_of(program)) %>% 
                                          mean()
                                  # Try to compute the test, it might fail if not enough cells are active
                                  data <- celltype_cell_score %>% dplyr::select(diagnosis, sample_name, phase, all_of(program)) %>% dplyr::filter(.[[program]]>0.1)
                                  colnames(data) <- c("diagnosis", "sample_name", "phase", "program")
                                  test_res <- try(lmer(program ~ diagnosis + phase + (1|sample_name), data = data),
                                                  silent = TRUE)
                                  if (inherits(test_res, "try-error")) {
                                          tests <- rbind(tests,
                                                         data_frame("cell_type"=cell_type, "program"=program, term = NA, "active_cells_control_perc"=perc_active_cells_control,
                                                                    "active_cells_case_perc"=perc_active_cells_case, "avg_control_activation"=avg_control_activation,
                                                                    "avg_case_activation"=avg_case_activation, "estimate"=NA, "std.err"=NA, "statistic"=NA, "p"=NA)
                                          )
                                  }else{
                                          test_res <- tidy(test_res, effects = "fixed")[2:3,]
                                          tests <- rbind(tests,
                                                         data_frame("cell_type"=cell_type, "program"=program, "term"=test_res$term,
                                                                    "active_cells_control_perc"=perc_active_cells_control, "active_cells_case_perc"=perc_active_cells_case,
                                                                    "avg_control_activation"=avg_control_activation, "avg_case_activation"=avg_case_activation,
                                                                    "estimate"=test_res$estimate, "std.err"=test_res$std.error, "statistic"=test_res$statistic, "p"=test_res$p.value)
                                          )
                                  }
                          }
                          return(tests)
                  },
                  mc.cores=snakemake@threads,
                  mc.preschedule=TRUE,
                  mc.cleanup=TRUE
)
GP_activation_lmm_phase_corrected <- bind_rows(res_lmm)
GP_activation_lmm_phase_corrected <- GP_activation_lmm_phase_corrected %>% 
        drop_na() %>% 
        group_by(cell_type) %>% 
        mutate(FDR = p.adjust(p, method = "fdr")) %>% 
        ungroup()
# Compute effect size as activationFC log2((avg_control_activation + estimate)/avg_case_activation)
GP_activation_lmm_phase_corrected <- GP_activation_lmm_phase_corrected %>% mutate(log2FC = log2((avg_control_activation+estimate)/avg_control_activation))

write_csv(GP_activation_lmm_phase_corrected, file = snakemake@output[["gp_activation_LMM"]])
message("Done")

# Close Logging ----------------------------------------------------------------
sink(type="message")
sink()