# ==============================================================================
# Script: rename_programs.R
# Author: Filippo Gastaldello
# Date: 018/08/2026
# Description: 
#   Tries to label gene programs detected by cNMF Performs ORA on the 
#   marker genes of the unknown gene programs.
#
# Snakemake Expected Inputs:
#   - snakemake@input[["cytopus_lsit"]] : Path to the cytopus gene sets json
#   - snakemake@input[["metadata]] : Path to seurat object's metadata
#   - Input paths for remaining inputs will be generated at runtime at the
#     start of the script
#
# Snakemake Expected Outputs:
#   - snakemake@output[["labeled_cell_scores"]] : Path for new cell score matrix
#                                                 with only labeled programs
#
# Snakemake Expected Params:
#   - snakemake@params[["output_dir"]] : Path of bse output directory
#   - snakemake@params[["k"]] : cNMF's k param, will be used to grep the needed
#                               files
#   - snakemake@params[["ldt"]] : cNMF's local density threshold, will be used 
#                                 to grep the needed files
#   - snakemake@params[["analysis_name"]] : Analysis name, will be used to grep
#                                           needed files
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

# 1. Generate paths ####
# ------------------------------------------------------------------------------
message("Generating input paths ...")
cell_scores_path <- paste0(snakemake@params[["output_dir"]],
                           "/",
                           snakemake@params[["analysis_name"]],
                           ".usages.k_",
                           snakemake@params[["k"]],
                           ".dt_",
                           gsub("\\.","_",snakemake@params[["ldt"]]),
                           ".consensus.txt"
                           )
markers_path <- paste0(snakemake@params[["output_dir"]],
                       "/",
                       snakemake@params[["analysis_name"]],
                       ".gene_spectra_score.k_",
                       snakemake@params[["k"]],
                       ".dt_",
                       gsub("\\.","_",snakemake@params[["ldt"]]),
                       ".txt")
message("Done")

# 2. Loading data ####
# ------------------------------------------------------------------------------
message("Loading cell score matrix from: ", cell_scores_path)
cell_scores <- read_tsv(cell_scores_path) %>% 
        dplyr::rename("cell_id"="...1")
message("Done")
message("Loading factor markers from: ", markers_path)
markers <- read_tsv(markers_path) %>% 
        dplyr::rename("program"="...1")
message("Done")
message("Loading cytopus list from : ", snakemake@input[["cytopus_list"]])
cytopus <- read_json(snakemake@input[["cytopus_list"]])
message("Done")
message("Loading seurat's metadata from: ", snakemake@input[["metadata"]])
metadata <- read_csv(snakemake@input[["metadata"]])
# Only keep metadata for cells in the spectra output
metadata <- metadata %>% dplyr::filter(cell_id %in% cell_scores$cell_id)
message("Done")

# 3. Annotate cell scores matrix with metadata ####
# ------------------------------------------------------------------------------
message("Annotating cell scores matrix with metadata ...")
# Add condition column to cell scores
cell_scores <- cell_scores %>% 
        left_join(metadata %>% dplyr::select(cell_id,
                                             diagnosis,
                                             sample_name,
                                             phase,
                                             celltype)) %>%
        relocate(diagnosis, .before = cell_id) %>%
        relocate(phase, .before = diagnosis) %>%
        relocate(sample_name, .before = diagnosis) %>% 
        relocate(celltype, .after=diagnosis)
message("Done")

# 4. Label unknown programs using ORA ####
# ------------------------------------------------------------------------------
message("Running Over Representation Analysis to label programs ...")
# Use Over Representation Analysis to give a name to programs

# Get name of programs to label
unknown_programs <- colnames(cell_scores)[5:ncol(cell_scores)]

# Refactor markers df
markers <- markers %>% pivot_longer(cols = 2:ncol(markers),
                                    names_to = "gene")

# Create background gene symbol universe
universe_symbol <- markers %>% 
        pull(gene) %>% 
        unique()

# Perform enrichment
res <- lapply(unknown_programs, function(factor){
        
        # Extract markers for current program
        current_markers <- markers %>% 
                                filter(program==factor) %>% 
                                arrange(desc(value)) %>% 
                                slice_head(n = 100) %>% 
                                pull(gene)
        
        cytopus_list <- unlist(cytopus, recursive = FALSE)
        names(cytopus_list) <- str_split_i(names(cytopus_list), "\\.",2)
        cytopus_list <- cytopus_list[!duplicated(names(cytopus_list))]
                             
        cytopus_term2gene <- data.frame(
                term = rep(names(cytopus_list), lengths(cytopus_list)),
                gene = unlist(cytopus_list, use.names = FALSE)
        )
        # Run ORA against Cytopus using factor markers
        ora_results <- enricher(current_markers,
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
                slice_min(p.adjust) %>% 
                pull(ID)
        # Keep track of which programs were named through ORA
        ora_results <- ifelse(is_empty(ora_results),
                              factor,
                              paste0(ora_results,"-(ORA)")) 
        return(data.frame(fact = factor,
                          name = paste0("cNMF-X-",
                                        factor,
                                        "-X-",
                                        ora_results)))
})
factor_names <- bind_rows(res)
# Add name to the programs that were identified
colnames(cell_scores)[6:ncol(cell_scores)] <- factor_names$name
# Identify programs that remained unknown and remove them
unknown_programs <- grep("\\d+$", colnames(cell_scores)[6:length(colnames(cell_scores))], value = TRUE)
cell_scores <- cell_scores %>%
        dplyr::select(-all_of(unknown_programs))
message("Done")

# 4. Save update cell scores matrix ####
# ------------------------------------------------------------------------------
message("Saving updated cell scores matrix to: ",
        snakemake@output[["labeled_cell_scores"]])
write_csv(cell_scores, snakemake@output[["labeled_cell_scores"]])
message("Done")

# Close Logging ----------------------------------------------------------------
sink(type="message")
sink()