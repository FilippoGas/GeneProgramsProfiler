# ==============================================================================
# Script: collect_benchmarks.R
# Author: Filippo Gastaldello
# Date: 03/11/2026
# Description:
#   Merge the per-rule Snakemake benchmarks of a single analysis into one
#   self-contained HTML table, rendered inline in the Snakemake report, and
#   summarize them with runtime and peak-memory bar charts.
#
# Snakemake Expected Inputs:
#   - snakemake@input[["benchmarks"]] : Benchmark .tsv files of every executed
#                                       rule, located at
#                                       results/<analysis>/<module>/benchmarks/
# Snakemake Expected Outputs:
#   - snakemake@output[["summary"]] : Destination path of the merged HTML table
#   - snakemake@output[["runtime"]] : Destination path of the runtime bar chart
#   - snakemake@output[["memory"]]  : Destination path of the memory bar chart
#
# ==============================================================================

# Setup Logging ----------------------------------------------------------------
# Redirect all output and messages to the Snakemake log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

# Prevent R from generating the default Rplots.pdf file
pdf(NULL)

# Load Libraries ---------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
})

# Helper Functions -------------------------------------------------------------
escape_html <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x
}

safe_median <- function(x) {
  if (all(is.na(x))) NA_real_ else stats::median(x, na.rm = TRUE)
}

numeric_metric <- function(bench, metric, i) {
  if (metric %in% names(bench)) {
    suppressWarnings(as.numeric(bench[[metric]][i]))
  } else {
    NA_real_
  }
}

empty_bar_plot <- function(x_label, title) {
  ggplot() +
    annotate("text", x = 0, y = 0, label = "No benchmark data available") +
    labs(x = x_label, y = NULL, title = title) +
    theme_bw(base_size = 11) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    )
}

bar_plot <- function(data, value, x_label, title) {
  data <- data[!is.na(data[[value]]) & data[[value]] > 0, , drop = FALSE]
  if (nrow(data) == 0) {
    return(empty_bar_plot(x_label, title))
  }
  ggplot(
    data,
    aes(
      x = .data[[value]],
      y = reorder(.data[["rule_label"]], .data[[value]]),
      fill = .data[["module"]]
    )
  ) +
    geom_col(width = 0.75) +
    scale_x_log10(labels = scales::label_number(accuracy = 0.1)) +
    labs(
      x = x_label,
      y = NULL,
      title = title,
      fill = "Module"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )
}

save_bar_plot <- function(plot, path) {
  n_bars <- length(unique(plot$data$rule_label))
  height <- max(4, 0.28 * n_bars + 1.6)
  ggsave(
    filename = path,
    plot = plot,
    device = "pdf",
    width = 9,
    height = height
  )
  message("Saved plot to ", path)
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

worker_pattern <- "^cNMF_factorize_worker_"

# Collect Benchmarks -----------------------------------------------------------
rows <- list()
tidy <- list()
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

    tidy[[length(tidy) + 1]] <- data.frame(
      module = module,
      rule = rule_name,
      s = numeric_metric(bench, "s", i),
      max_rss = numeric_metric(bench, "max_rss", i),
      stringsAsFactors = FALSE
    )
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

# Summarize Benchmarks ---------------------------------------------------------
# Collapse the many cNMF factorization workers into a single median row so the
# charts stay readable; all other rules are already one row per file.
bench <- if (length(tidy) > 0) {
  dplyr::bind_rows(tidy)
} else {
  data.frame(
    module = character(0),
    rule = character(0),
    s = numeric(0),
    max_rss = numeric(0)
  )
}

bench$rule_label <- as.character(ifelse(
  grepl(worker_pattern, bench$rule),
  "cNMF_factorize_worker",
  bench$rule
))

plot_data <- bench %>%
  dplyr::group_by(.data[["module"]], .data[["rule_label"]]) %>%
  dplyr::summarise(
    s = safe_median(.data[["s"]]),
    max_rss = safe_median(.data[["max_rss"]]),
    n = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    rule_label = dplyr::if_else(
      .data[["n"]] > 1,
      paste0(.data[["rule_label"]], " (median of ", .data[["n"]], ")"),
      .data[["rule_label"]]
    )
  )

runtime_plot <- bar_plot(
  plot_data,
  value = "s",
  x_label = "Wall-clock time [s]",
  title = "Runtime per rule"
)
save_bar_plot(runtime_plot, snakemake@output[["runtime"]])

memory_plot <- bar_plot(
  plot_data,
  value = "max_rss",
  x_label = "Peak resident memory [MiB]",
  title = "Peak memory per rule"
)
save_bar_plot(memory_plot, snakemake@output[["memory"]])

message("Done")
