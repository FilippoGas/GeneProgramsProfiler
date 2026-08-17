# ==============================================================================
# Script: extract_best_k.R
# Author: Filippo Gastaldello
# Date: 14/08/2026
# Description: 
#   Find the best value of k (gene programs to look for in the dataset) given 
#   the results of the previous rules
#
# Snakemake Expected Inputs:
#   - snakemake@input[["k_selection_stats"]] : Path to dataframe holding stats
#                                              about k selection
# Snakemake Expected Outputs:
#   - snakemake@output[["best_k"]] : Path to save the best k to a file 
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
        library(reticulate)
})
# Import numpy
np <- import("numpy")

message("Starting R script \"extract_best_k.R\"...")

# 1. Load Data ####
# ------------------------------------------------------------------------------
message("Loading stats on k selection from: ",
        snakemake@input[["k_selection_stats"]])
k_stats <- np$load(snakemake@input[["k_selection_stats"]], allow_pickle=TRUE)
message("Done")

# 2. Choose best K ####
# ------------------------------------------------------------------------------
message("Extracting data from compressed numpy dataframe ...")
df_data   <- k_stats$f[["data"]]
df_index  <- k_stats$f[["index"]]
df_column <- k_stats$f[["columns"]]
# 2. Convert the 2D data array into an R data frame
my_dataframe <- as.data.frame(df_data)
# 3. Assign the row and column names
# We wrap them in as.character() to ensure they are properly formatted as text vectors
rownames(my_dataframe) <- as.character(df_index)
colnames(my_dataframe) <- as.character(df_column)
message("Done")
message("Computing best K ...")
# Scale stability and error in [0.1] and compute delta
my_dataframe <- my_dataframe %>%
        mutate(silhouette=(silhouette-min(silhouette))/(max(silhouette)-min(silhouette)),
               prediction_error=(prediction_error-min(prediction_error))/(max(prediction_error)-min(prediction_error)),
               delta=silhouette-prediction_error)
# Always close the connection
k_stats$close()
message("Done")
message("Saving best K to a :", snakemake@output["best_k"])
selected_k <- my_dataframe[my_dataframe$delta==max(my_dataframe$delta),"k"]
write(as.character(selected_k), snakemake@output["best_k"])
message("Done")

# Close Logging ----------------------------------------------------------------
sink(type="message")
sink()