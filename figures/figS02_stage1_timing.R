#!/usr/bin/env Rscript
# Fig S2: matched Stage 1 cost: fit + SE time, whole-process time and peak memory by workload and
# method (frozen S2 timing means). Reads results/P17_simulations.
suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(data.table)})
.repo <- system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE)
source(file.path(.repo, "figures", "theme.R"))
tm <- as.data.table(read_result("08_simulations/speed_benchmark", "stage1_panel_times"))
ds <- as.data.table(read_result("08_simulations/summary", "table_s1_design_provenance"))[study == "S2"]
mlab <- c(ordinary_rif_ols = "RIF-OLS", gmoq = "G-MoQ, pooled", gmoq_jackknife = "G-MoQ, jackknife")
mcol <- setNames(c(paper_palette[["regenie"]], paper_palette[["gmoq"]], "#7fb3d5"), mlab)
tm[, method_f := factor(mlab[method], levels = mlab)]
wl <- unique(tm[, .(workload, n, m)])[order(n, m)]
fmt <- function(v) format(v, big.mark = ",", scientific = FALSE, trim = TRUE)
wl[, lab := sprintf("%s\n%s", fmt(n), fmt(m))]
tm[, wl_f := factor(wl$lab[match(workload, wl$workload)], levels = wl$lab)]
mtr <- c("fit + SE seconds", "whole-process seconds", "peak RSS MiB")
ylab <- c("Fit + SE time (s)", "Whole-process time (s)", "Peak memory (MiB)")
# One panel per metric (a-c); the metric is the y-axis title. The key sits inside panel a.
metric_panel <- function(i) {
  ggplot(tm[metric == mtr[i]], aes(wl_f, mean_value, colour = method_f)) +
    geom_point(position = position_dodge(width = 0.55), size = 1.6) +
    scale_y_log10(labels = function(v) format(v, big.mark = ",", scientific = FALSE, drop0trailing = TRUE, trim = TRUE)) +
    scale_colour_manual(values = mcol, name = NULL, drop = FALSE, guide = if (i == 1) "legend" else "none") +
    labs(x = if (i == 2) "Sample size and variants per run" else NULL, y = ylab[i], tag = letters[i]) +
    theme_paper() + legend_inside(0.02, 0.98, c(0, 1))
}
p <- wrap_plots(lapply(1:3, metric_panel), nrow = 1)
reps <- unique(tm$count)
cap <- paste(
  "Figure S2. Computational cost of Stage 1 G-MoQ and ordinary RIF-OLS on matched workloads.",
  sprintf(paste("Four simulated workloads, each a sample size (top tick label) and a number of variants fitted in one run",
                "(bottom tick label): %s. Each run fits all 45 quantile levels on %d core%s."),
          paste(wl[, .(s = sprintf("%s participants with %s variants", fmt(n[1]), paste(fmt(m), collapse = " or "))), by = n]$s,
                collapse = ", and "),
          as.integer(ds$threads), if (as.integer(ds$threads) == 1) "" else "s"),
  "Methods are ordinary RIF-OLS fitted separately at each quantile level, and G-MoQ Stage 1 with the pooled",
  "covariance or with the 1,000-block jackknife covariance (Methods).",
  sprintf(paste("Each point is the arithmetic mean of %s timing repeats of one method on one workload, on a log scale;",
                "no error bars are shown. **a**, Time to fit the coefficients and their standard errors (s).",
                "**b**, Wall time of the whole process (s). **c**, Peak resident memory (MiB)."),
          if (length(reps) == 1 && reps == 3) "three" else paste(reps, collapse = "/")),
  "No statistical tests are shown.",
  "RIF-OLS, recentred influence function ordinary least squares; SE, standard error; MiB, mebibyte.")
save_figure(p, "figS02_stage1_timing", tm[, .(workload, n, m, method, metric, mean_value, timing_repeats = count)], cap,
            width = 183, height = 65,
            results_used = c("08_simulations/speed_benchmark", "08_simulations/summary"))
