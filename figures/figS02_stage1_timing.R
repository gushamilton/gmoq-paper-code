#!/usr/bin/env Rscript
# Supplementary Figure S2: computational cost of Stage 1 against ordinary RIF-OLS
#
# What: on four matched simulated workloads (a sample size and a number of variants per run),
#   compares the time and memory of ordinary RIF-OLS fitted separately at each quantile level
#   with G-MoQ Stage 1 (pooled or jackknife covariance). Each point is a mean over timing
#   repeats.
#   a: time to fit coefficients and standard errors.
#   b: wall time of the whole process.
#   c: peak resident memory.
# Paper: Supplementary Fig. S2
# Input: results/08_simulations/summary/{s2_mean_timings,table_s1_design_provenance}.csv
# Output: paper_assets/figS02_stage1_timing.{png,pdf}, _source.csv, _caption.md
# Runs: locally from results/ in about 2 seconds

suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(data.table)})
source(file.path(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE), "figures", "theme.R"))

## 1. Load ----------------------------------------------------------------------------------
timings <- as.data.table(read_result("08_simulations/summary", "s2_mean_timings"))
design <- as.data.table(read_result("08_simulations/summary", "table_s1_design_provenance"))[study == "S2"]
method_lab <- c(ordinary_rif_ols = "RIF-OLS", gmoq = "G-MoQ, pooled", gmoq_jackknife = "G-MoQ, jackknife")
method_col <- setNames(c(paper_palette[["regenie"]], paper_palette[["gmoq"]], "#7fb3d5"), method_lab)
timings[, method_f := factor(method_lab[method], levels = method_lab)]
fmt <- function(v) format(v, big.mark = ",", scientific = FALSE, trim = TRUE)
# Workload tick labels: sample size over number of variants.
workloads <- unique(timings[, .(workload, n, m)])[order(n, m)]
workloads[, lab := sprintf("%s\n%s", fmt(n), fmt(m))]
timings[, wl_f := factor(workloads$lab[match(workload, workloads$workload)], levels = workloads$lab)]

## 2. Panels a-c: one panel per cost measure ----------------------------------------------
metrics <- c("fit + SE seconds", "whole-process seconds", "peak RSS MiB")
y_titles <- c("Fit + SE time (s)", "Whole-process time (s)", "Peak memory (MiB)")
metric_panel <- function(i) {
  ggplot(timings[metric == metrics[i]], aes(wl_f, mean_value, colour = method_f)) +
    geom_point(position = position_dodge(width = 0.55), size = 1.6) +
    scale_y_log10(labels = function(v) format(v, big.mark = ",", scientific = FALSE, drop0trailing = TRUE, trim = TRUE)) +
    scale_colour_manual(values = method_col, name = NULL, drop = FALSE, guide = if (i == 1) "legend" else "none") +
    labs(x = if (i == 2) "Sample size and variants per run" else NULL, y = y_titles[i], tag = letters[i]) +
    theme_paper() + legend_inside(0.02, 0.98, c(0, 1))
}

## 3. Assemble and save ---------------------------------------------------------------------
fig <- wrap_plots(lapply(1:3, metric_panel), nrow = 1)
reps <- unique(timings$count)
threads <- as.integer(design$threads)
cap <- paste(
  "Figure S2. Computational cost of Stage 1 G-MoQ and ordinary RIF-OLS on matched workloads.",
  sprintf(paste("Four simulated workloads, each a sample size (top tick label) and a number of variants fitted in one run",
                "(bottom tick label): %s. Each run fits all 45 quantile levels on %d core%s."),
          paste(workloads[, .(s = sprintf("%s participants with %s variants", fmt(n[1]), paste(fmt(m), collapse = " or "))), by = n]$s,
                collapse = ", and "),
          threads, if (threads == 1) "" else "s"),
  "Methods are ordinary RIF-OLS fitted separately at each quantile level, and G-MoQ Stage 1 with the pooled",
  "covariance or with the 1,000-block jackknife covariance (Methods).",
  sprintf(paste("Each point is the arithmetic mean of %s timing repeats of one method on one workload, on a log scale;",
                "no error bars are shown. **a**, Time to fit the coefficients and their standard errors (s).",
                "**b**, Wall time of the whole process (s). **c**, Peak resident memory (MiB)."),
          if (length(reps) == 1 && reps == 3) "three" else paste(reps, collapse = "/")),
  "No statistical tests are shown.",
  "RIF-OLS, recentred influence function ordinary least squares; SE, standard error; MiB, mebibyte.")
save_figure(fig, "figS02_stage1_timing", timings[, .(workload, n, m, method, metric, mean_value, timing_repeats = count)], cap,
            width = 183, height = 65, results_used = "08_simulations/summary")
