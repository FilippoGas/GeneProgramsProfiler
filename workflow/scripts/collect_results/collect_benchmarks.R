# ==============================================================================
# Script: collect_benchmarks.R
# Author: Filippo Gastaldello
# Date: 03/11/2026
# Description:
#   Merge the per-rule Snakemake benchmarks of a single analysis into one
#   self-contained HTML table, rendered inline in the Snakemake report.
#
# Snakemake Expected Inputs:
#   - snakemake@input[["benchmarks"]] : Benchmark .tsv files of every executed
#                                       rule, located at
#                                       results/<analysis>/<module>/benchmarks/
# Snakemake Expected Outputs:
#   - snakemake@output[["summary"]] : Destination path of the merged HTML table
#
# ==============================================================================

# Setup Logging ----------------------------------------------------------------
# Redirect all output and messages to the Snakemake log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

# Helper Functions -------------------------------------------------------------
escape_html <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x
}

# Constants --------------------------------------------------------------------
metrics <- c(
  "s", "h:m:s", "max_rss", "max_vms", "max_uss",
  "max_pss", "io_in", "io_out", "mean_load", "cpu_time"
)

style <- paste0(
  "body{font-family:sans-serif;margin:0;padding:1rem}",
  "table{border-collapse:collapse;font-size:.85rem}",
  "th,td{border:1px solid #cbd5e1;padding:.25rem .5rem;text-align:right}",
  "th{background:#1e293b;color:#f8fafc}",
  "th:nth-child(-n+2),td:nth-child(-n+2){text-align:left}"
)

# Collect Benchmarks -----------------------------------------------------------
rows <- list()
for (path in snakemake@input[["benchmarks"]]) {
  module <- basename(dirname(dirname(path)))
  rule_name <- tools::file_path_sans_ext(basename(path))
  bench <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)

  for (i in seq_len(nrow(bench))) {
    values <- vapply(
      metrics,
      function(metric) {
        if (metric %in% names(bench)) as.character(bench[[metric]][i]) else ""
      },
      character(1)
    )
    rows[[length(rows) + 1]] <- c(module, rule_name, values)
  }
}

if (length(rows) == 0) {
  rows <- list(c("", "", rep("", length(metrics))))
}

modules <- vapply(rows, function(row) row[1], character(1))
rules <- vapply(rows, function(row) row[2], character(1))
rows <- rows[order(tolower(modules), rules)]

# Render HTML ------------------------------------------------------------------
header <- paste0("<th>", escape_html(metrics), "</th>", collapse = "")
body <- vapply(
  rows,
  function(row) {
    paste0(
      "<tr>",
      paste0("<td>", escape_html(row), "</td>", collapse = ""),
      "</tr>"
    )
  },
  character(1)
)

html <- c(
  "<!DOCTYPE html>",
  "<html>",
  "<head>",
  "<meta charset='utf-8'>",
  paste0("<style>", style, "</style>"),
  "</head>",
  "<body>",
  "<table>",
  paste0("<thead><tr><th>module</th><th>rule</th>", header, "</tr></thead>"),
  "<tbody>",
  body,
  "</tbody>",
  "</table>",
  "</body>",
  "</html>"
)

writeLines(html, snakemake@output[["summary"]])

message("Merged ", length(rows), " benchmark rows into ", snakemake@output[["summary"]])
