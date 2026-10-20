#!/usr/bin/env Rscript
# Supplementary Figure S5: computational cost of G-MoQ and direct mixture maximum likelihood
#
# What: compares the cost of one matched panel of simulated variants under the two-component
#   mixture model for the whole G-MoQ pipeline (including one-off costs) and for refitting the
#   mixture by maximum likelihood at every variant, then scales the per-variant costs of
#   Fig. 2c to a genome-wide analysis.
#   a: time to fit estimates and standard errors.
#   b: wall time of the whole process.
#   c: peak resident memory.
#   d: implied CPU-hours for 1.2 million variants at 400,000 participants.
# Paper: Supplementary Fig. S5
# Input: results/08_simulations/summary/s5_mean_timings.csv,
#   results/03_testosterone/speed/per_snp_costs.csv
# Output: paper_assets/figS05_mixture_timing.{png,pdf}, _source.csv, _caption.md
# Runs: locally from results/ in about 2 seconds

suppressPackageStartupMessages({library(ggplot2); library(data.table); library(patchwork)})
source(file.path(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE), "figures", "theme.R"))

fmt <- function(v) format(v, big.mark = ",", scientific = FALSE, trim = TRUE)

## 1. Panels a-c: one matched panel, one panel per cost measure -----------------------------
timings <- as.data.table(read_result("08_simulations/summary", "s5_mean_timings"))
method_lab <- c(gmoq_pooled = "G-MoQ, pooled", gmoq_jackknife = "G-MoQ, jackknife", mixture_ml = "Mixture ML")
method_col <- setNames(c(paper_palette[["gmoq"]], "#7fb3d5", paper_palette[["mixture"]]), method_lab)
timings[, method_f := factor(method_lab[method], levels = method_lab)]
metrics <- c("fit + SE seconds", "whole-process seconds", "peak RSS MiB")
y_titles <- c("Fit + SE time (s)", "Whole-process time (s)", "Peak memory (MiB)")
metric_panel <- function(i) {
  ggplot(timings[metric == metrics[i]], aes(n, mean_value, colour = method_f)) +
    geom_line(linewidth = 0.7) + geom_point(size = 1.4) +
    scale_x_log10(breaks = c(1e4, 5e4, 4e5), labels = c("10,000", "50,000", "400,000")) +
    scale_y_log10(labels = function(v) format(v, big.mark = ",", scientific = FALSE, drop0trailing = TRUE, trim = TRUE)) +
    scale_colour_manual(values = method_col, name = NULL, drop = FALSE, guide = if (i == 1) "legend" else "none") +
    labs(x = "Sample size", y = y_titles[i], tag = letters[i]) +
    theme_paper() + legend_inside(0.02, 0.98, c(0, 1))
}
pABC <- wrap_plots(lapply(1:3, metric_panel), nrow = 1)

## 2. Panel d: genome-wide CPU-hours at N = 400,000 ------------------------------------------
# Per-variant costs times 1.2 million variants; bars start from a common floor of 0.01 h.
costs <- as.data.table(read_result("03_testosterone/speed", "per_snp_costs"))[n == 4e5]
bar_lab <- c("Direct mixture ML", "G-MoQ Stage 1 + Stage 2", "Stage 2 alone")
bar_col <- setNames(c(paper_palette[["mixture"]], paper_palette[["gmoq"]], "#2a8c82"), bar_lab)
n_variants <- 1.2e6
gw <- data.table(what = factor(bar_lab, bar_lab),
                 per_snp = c(costs$ml_per_snp, costs$stage1_per_snp + costs$stage2_per_snp_param_gwas,
                             costs$stage2_per_snp_param_gwas))
gw[, cpu_h := per_snp * n_variants / 3600]
gw[, lab := ifelse(cpu_h >= 100, format(round(cpu_h), big.mark = ","), ifelse(cpu_h >= 1, sprintf("%.1f", cpu_h), sprintf("%.2f", cpu_h)))]
pD <- ggplot(gw, aes(y = what, colour = what)) +
  geom_segment(aes(x = 0.01, xend = cpu_h, yend = what), linewidth = 6, lineend = "butt") +
  geom_text(aes(x = cpu_h, label = lab), hjust = -0.15, size = 6 / .pt, colour = "grey15") +
  scale_x_log10(breaks = c(0.01, 0.1, 1, 10, 100, 1000), labels = c("0.01", "0.1", "1", "10", "100", "1,000"),
                limits = c(0.01, NA), expand = expansion(mult = c(0, 0.12))) +
  scale_y_discrete(limits = rev(bar_lab)) + scale_colour_manual(values = bar_col, guide = "none") +
  labs(x = "CPU-hours on one core (log scale)", y = NULL, tag = "d") +
  theme_paper() + theme(panel.grid.major.y = element_blank(), plot.margin = margin(4, 12, 4, 4))

## 3. Assemble and save ---------------------------------------------------------------------
fig <- pABC / free(pD) + plot_layout(heights = c(1.4, 1))
reps <- unique(timings$count)
cap <- paste(
  "Figure S5. Computational cost of G-MoQ and direct mixture maximum likelihood under the two-component mixture model.",
  sprintf(paste("Panels a-c give the cost of one matched panel of %s simulated variants at three sample sizes, for the",
                "whole G-MoQ pipeline with the pooled or the 1,000-block jackknife covariance and for direct mixture",
                "maximum likelihood. **a**, Time to fit the estimates and their standard errors (s). **b**, Wall time of the",
                "whole process (s). **c**, Peak resident memory (MiB). Each G-MoQ time includes the one-off baseline mixture fit,",
                "Stage 1 and Stage 2 for the panel, so it is dominated by costs not repeated per variant (per-variant",
                "costs, Fig. 2c). Each point is the arithmetic mean of %s timing repeats of one method at one sample size,",
                "on log axes; no error bars are shown."),
          unique(timings$M), if (length(reps) == 1 && reps == 3) "three" else paste(reps, collapse = "/")),
  sprintf(paste("**d**, Implied cost of %s million variants at %s participants on one core, from the per-variant",
                "costs of Fig. 2c; the value printed is CPU-hours, and bars start at 0.01. Direct mixture maximum likelihood",
                "needs %s CPU-hours; G-MoQ Stages 1 and 2, %.1f CPU-hours (%.0f times faster); Stage 2 alone given",
                "the Stage 1 output, %.0f CPU-minutes (about %s times faster). One-off per-trait costs are excluded."),
          n_variants / 1e6, fmt(4e5), fmt(round(gw$cpu_h[1])), gw$cpu_h[2], gw$cpu_h[1] / gw$cpu_h[2], gw$cpu_h[3] * 60,
          fmt(round(gw$cpu_h[1] / gw$cpu_h[3], -3))),
  "No statistical tests are shown.",
  "ML, maximum likelihood; SE, standard error; MiB, mebibyte; CPU, central processing unit.")
src <- rbind(timings[, .(panel = "abc", n, M, method, metric, value = mean_value, timing_repeats = count)],
             gw[, .(panel = "d", n = 4e5, M = n_variants, method = as.character(what), metric = "genome-wide CPU-hours",
                    value = cpu_h, timing_repeats = NA_integer_)])
save_figure(fig, "figS05_mixture_timing", src, cap, width = 183, height = 110,
            results_used = c("08_simulations/summary", "03_testosterone/speed"))
