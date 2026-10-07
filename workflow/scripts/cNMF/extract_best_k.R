# ==============================================================================
# Script: extract_best_k.R
# Author: Filippo Gastaldello
# Date: 14/08/2026
# Description:
#   Find the best value of k (gene programs to look for in the dataset) given
#   the results of the previous rules. Rows whose k-selection statistics are
#   broken (non-finite or exactly-zero prediction_error) are discarded rather
#   than allowed to win the ranking.
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
# source() buffers warnings until it returns, i.e. after the sinks below are
# closed; print them as they occur so they end up in the Snakemake log instead
# of on stderr where nobody looks.
options(warn = 1)

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
# Always close the connection
k_stats$close()
message("Done")

# 2. Validate the raw statistics ####
# ------------------------------------------------------------------------------
# prediction_error is the squared residual between the refit usages and the
# normalized counts, so a value of exactly 0 (or any non-finite value) can never
# mark a genuinely good k: it means the k-selection run failed for that k.
# delta is built from a min-max scaling of this column, and a bogus 0 would
# scale to the *best* possible error and win the ranking, silently forcing the
# downstream consensus onto an absurd k. Discard such rows instead.
usable <- is.finite(my_dataframe$prediction_error) &
        my_dataframe$prediction_error > 0 &
        is.finite(my_dataframe$silhouette)

if (!any(usable)) {
        stop("No k has usable statistics (prediction_error must be finite and ",
             "> 0, silhouette must be finite) in ",
             snakemake@input[["k_selection_stats"]],
             ". Re-run cNMF_k_selection_plot; k selection is impossible.",
             call. = FALSE)
}
if (any(!usable)) {
        warning("Discarding k = ",
                paste(my_dataframe$k[!usable], collapse = ", "),
                ": prediction_error is non-finite or exactly 0, or ",
                "silhouette is non-finite. These k cannot be selected. ",
                "This indicates cNMF's k_selection_plot failed for them.",
                call. = FALSE)
}
if (sum(usable) == 1) {
        warning("Only k = ", my_dataframe$k[usable],
                " has usable statistics; k selection is degenerate.",
                call. = FALSE)
}

# 3. Compute delta ####
# ------------------------------------------------------------------------------
message("Computing best K ...")
# Scale stability and error in [0,1] and compute delta.
# Only usable rows define each column's range, so a discarded row can never
# shift the ranking of the good ones. A zero-width range means the column
# carries no discriminating information and is scaled to a constant 0.
scale01 <- function(x, valid) {
        rng <- range(x[valid])
        if (!all(is.finite(rng)) || diff(rng) == 0) {
                return(rep(0, length(x)))
        }
        (x - rng[1]) / diff(rng)
}
my_dataframe$delta <- scale01(my_dataframe$silhouette, usable) -
        scale01(my_dataframe$prediction_error, usable)
my_dataframe$delta[!usable] <- -Inf
message("Done")

# 4. Save best K ####
# ------------------------------------------------------------------------------
# which.max() returns the first maximum, so ties are broken deterministically
# toward the smallest k (the more conservative choice) because cNMF writes the
# stats sorted by ascending k.
best_idx <- which.max(my_dataframe$delta)
selected_k <- my_dataframe$k[best_idx]
n_tied <- sum(my_dataframe$delta == my_dataframe$delta[best_idx])
if (n_tied > 1) {
        message("Tie in delta between k = ",
                paste(my_dataframe$k[my_dataframe$delta == my_dataframe$delta[best_idx]],
                      collapse = ", "),
                "; selecting the smallest k = ", selected_k)
}
message("Saving best K to :", snakemake@output[["best_k"]])
write(x = as.character(selected_k), file = snakemake@output[["best_k"]])
message("Done")

# Close Logging ----------------------------------------------------------------
sink(type="message")
sink()
