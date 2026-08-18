# ==============================================================================
# Script: rename_programs.R
# Author: Filippo Gastaldello
# Date: 06/08/2026
# Description: 
#   Tries to label factor there weren't labeled by spectra. Performs ORA on the 
#   marker genes of the unknown factors.
#
# Snakemake Expected Inputs:
#   - snakemake@input[["cell_scores"]] : Path to the cell score matrix
#                                        (spectra's output)
#   - snakemake@input[["factor_markers"]] : Path to the factors markers matrix
#                                          (spectra's output)
#   - snakemake@input[["cytopus_lsit"]] : Path to the cytopus gene sets json
#   - snakemake@input[["metadata]] : Path to seurat object's metadata
#
# Snakemake Expected Outputs:
#   - snakemake@output[["cell_scores"]] : Path for new cell score matrix, with
#                                         only labeled programs
#                                       
# ==============================================================================

# Setup Logging ----------------------------------------------------------------
# Redirect all output and messages to the Snakemake log file
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# Load Libraries ---------------------------------------------------------------
suppressPackageStartupMessages({
        library(tidyverse)
        library(jsonlite)
        library(clusterProfiler)
        library(org.Hs.eg.db)
})

message("Starting R script \"rename_programs.R\"...")

# 1. Load Data ####
# ------------------------------------------------------------------------------
message("Loading cell score matrix from: ", snakemake@input[["cell_scores"]])
cell_scores <- read_csv(snakemake@input[["cell_scores"]]) %>% 
        dplyr::rename("cell_id"="...1")
message("Done")
message("Loading factor markers from: ", snakemake@input[["factor_markers"]])
markers <- read_csv(snakemake@input[["factor_markers"]])
message("Done")
message("Loading cytopus list from : ", snakemake@input[["cytopus_list"]])
cytopus <- read_json(snakemake@input[["cytopus_list"]])
message("Done")
message("Loading seurat's metadata from: ", snakemake@input[["metadata"]])
metadata <- read_csv(snakemake@input[["metadata"]])
# Only keep metadata for cells in the spectra output
metadata <- metadata %>% dplyr::filter(cell_id %in% cell_scores$cell_id)
message("Done")

# 2. Annotate cell scores matrix with metadata ####
# ------------------------------------------------------------------------------
message("Annotating cell scores matrix with metadata ...")
# Add condition column to cell scores
cell_scores <- cell_scores %>% 
        left_join(metadata %>% dplyr::select(cell_id,
                                             diagnosis,
                                             phase,
                                             sample_name,
                                             celltype)) %>%
        relocate(diagnosis, .before = cell_id) %>%
        relocate(phase, .after = diagnosis) %>%
        relocate(sample_name, .before = diagnosis) %>% 
        relocate(celltype, .after=diagnosis)
message("Done")

# 3. Label unknown programs using ORA ####
# ------------------------------------------------------------------------------
message("Splitting cell scores matrix in known and unknown factors ...")
# Split cell_score in known and unknown programs
unknown_programs <- grep("\\d+$",
                         colnames(cell_scores)[6:length(colnames(cell_scores))],
                         value = TRUE)
unknown_cell_scores <- cell_scores %>%
        dplyr::select(all_of(c(colnames(cell_scores)[1:5],
                               unknown_programs)
                             )
                      )
cell_scores <- cell_scores %>%
        dplyr::select(-all_of(unknown_programs))
message("Done")

message("Running Over Representation Analysis to label unknown factors ...")
# Use Over Representation Analysis to give a name to unknown programs
# Get factor markers to use in the enrichment 
markers <- markers %>% dplyr::select(all_of(unknown_programs))

# Create background gene symbol universe
universe_symbol <- markers %>% 
        pivot_longer(everything()) %>% 
        pull(value) %>% 
        unique()

# Perform enrichment
res <- lapply(unknown_programs, function(factor){
        
        cell_type <- str_split_i(factor, "-X-", 2)
        cell_type <- ifelse(cell_type == "global", "all", cell_type)
        
        # Prepare Cytopus term2gene for current cell type (program specificity)
        if(cell_type=="all") {
                cytopus_list <- c(cytopus[["global"]])
        }else{
                cytopus_list <- c(cytopus[[cell_type]],
                                  cytopus[["global"]])
        }
        cytopus_list <- lapply(cytopus_list,
                               function(x) as.character(unlist(x,
                                                               use.names = FALSE
                                                               )
                                                        )
                               )
        cytopus_term2gene <- data.frame(
                term = rep(names(cytopus_list), lengths(cytopus_list)),
                gene = unlist(cytopus_list, use.names = FALSE)
        )
        # Run ORA against Cytopus using factor markers
        ora_results <- enricher(markers %>% pull(factor),
                                minGSSize = 5,
                                universe  = universe_symbol,
                                TERM2GENE = cytopus_term2gene)
        
        if(is.null(ora_results)){
                return(data.frame(fact = factor,
                                  name = factor)
                       )
                }
        
        ora_results <- ora_results@result %>% 
                filter(p.adjust <= 0.05) %>% 
                mutate(GeneRatio = as.numeric(str_split_i(GeneRatio, "/",1))/
                               as.numeric(str_split_i(GeneRatio, "/",2))) %>% 
                slice_min(GeneRatio) %>% 
                pull(ID)
        # Keep track of which programs were named through ORA
        ora_results <- ifelse(is_empty(ora_results),
                              str_split_i(factor, "-X-",1),
                              paste0(ora_results,"-(ORA)")) 
        return(data.frame(fact = factor,
                          name = paste0(str_split_i(factor, "-X-",1),
                                        "-X-",
                                        str_split_i(factor, "-X-",2),
                                        "-X-",
                                        ora_results)))
})
factor_names <- bind_rows(res)
# Add name to the programs that were identified
colnames(unknown_cell_scores)[6:ncol(unknown_cell_scores)] <- factor_names$name
# Identify programs that remained unknown and remove them
unknown_programs <- grep("\\d+$", colnames(unknown_cell_scores)[6:length(colnames(unknown_cell_scores))], value = TRUE)
unknown_cell_scores <- unknown_cell_scores %>%
        dplyr::select(-all_of(unknown_programs))
# Merge cell_scores and unknown_cell_scores
message(colnames(cell_scores))
message(colnames(unknown_cell_scores))
cell_scores <- cell_scores %>%
        left_join(unknown_cell_scores %>% dplyr::select(-c(1:3,5)),
                  by = "cell_id")
message("Done")

# 4. Save update cell scores matrix ####
# ------------------------------------------------------------------------------
message("Saving updated cell scores matrix to: ",
        snakemake@output[["cell_scores"]])
write_csv(cell_scores, snakemake@output[["cell_scores"]])
message("Done")

# Close Logging ----------------------------------------------------------------
sink(type="message")
sink()