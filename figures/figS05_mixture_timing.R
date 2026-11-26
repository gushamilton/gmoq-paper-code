#!/usr/bin/env Rscript
# Fig S5: complete matched M = 100 mixture panel cost (fit + SE, whole-process time, peak memory)
# for the whole G-MoQ pipeline (pooled, jackknife) and full mixture ML by N (frozen S5 timing means);
# (d) genome-wide CPU-hours at N = 400,000 from the per-SNP costs of aux_speed_profile (as Fig 2C).
suppressPackageStartupMessages({library(ggplot2); library(data.table); library(patchwork)})
.repo <- system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE)
source(file.path(.repo, "figures", "theme.R"))
tm_w <- as.data.table(read_result("08_simulations/speed_benchmark", "mixture_panel_times"))[method %in% c("gmoq_pooled", "gmoq_jackknife", "mixture_ml")]
tm <- rbind(
  tm_w[, .(n, M = variants_per_panel, method, metric = "fit + SE seconds",      mean_value = mean_fit_seconds,  count = panels)],
  tm_w[, .(n, M = variants_per_panel, method, metric = "whole-process seconds", mean_value = mean_wall_seconds, count = panels)],
  tm_w[, .(n, M = variants_per_panel, method, metric = "peak RSS MiB",          mean_value = mean_peak_rss_mib, count = panels)]
)
mlab <- c(gmoq_pooled = "G-MoQ, pooled", gmoq_jackknife = "G-MoQ, jackknife", mixture_ml = "Mixture ML")
mcol <- setNames(c(paper_palette[["gmoq"]], "#7fb3d5", paper_palette[["mixture"]]), mlab)
tm[, method_f := factor(mlab[method], levels = mlab)]
mtr <- c("fit + SE seconds", "whole-process seconds", "peak RSS MiB")
ylab <- c("Fit + SE time (s)", "Whole-process time (s)", "Peak memory (MiB)")
# One panel per metric (a-c); the metric is the y-axis title. The key sits inside panel a.
metric_panel <- function(i) {
  ggplot(tm[metric == mtr[i]], aes(n, mean_value, colour = method_f)) +
    geom_line(linewidth = 0.7) + geom_point(size = 1.4) +
    scale_x_log10(breaks = c(1e4, 5e4, 4e5), labels = c("10,000", "50,000", "400,000")) +
    scale_y_log10(labels = function(v) format(v, big.mark = ",", scientific = FALSE, drop0trailing = TRUE, trim = TRUE)) +
    scale_colour_manual(values = mcol, name = NULL, drop = FALSE, guide = if (i == 1) "legend" else "none") +
    labs(x = "Sample size", y = ylab[i], tag = letters[i]) +
    theme_paper() + legend_inside(0.02, 0.98, c(0, 1))
}
pABC <- wrap_plots(lapply(1:3, metric_panel), nrow = 1)

# ---- d: projected genome-wide CPU-hours at N = 400,000 (per-SNP costs x 1.2 million SNPs)
pc <- as.data.table(read_result("08_simulations/speed_benchmark", "component_costs"))[n == 4e5]
d_lv <- c("Direct mixture ML", "G-MoQ Stage 1 + Stage 2", "Stage 2 alone")
d_col <- setNames(c(paper_palette[["mixture"]], paper_palette[["gmoq"]], "#2a8c82"), d_lv)  # teal, distinct from the jackknife light blue
G_snps <- 1.2e6
g <- data.table(what = factor(d_lv, d_lv),
                per_snp = c(pc$direct_seconds_per_variant,
                            pc$stage1_seconds_per_extra_variant + pc$stage2_seconds_per_variant_given_summaries,
                            pc$stage2_seconds_per_variant_given_summaries))
g[, cpu_h := per_snp * G_snps / 3600]
g[, lab := ifelse(cpu_h >= 100, sprintf("%s h", format(round(cpu_h), big.mark = ",")),
           ifelse(cpu_h >= 1, sprintf("%.2f h", cpu_h), sprintf("%.1f s", cpu_h * 3600)))]
pD <- ggplot(g, aes(y = what, colour = what)) +
  geom_segment(aes(x = 0.001, xend = cpu_h, yend = what), linewidth = 6, lineend = "butt") +
  geom_text(aes(x = cpu_h, label = lab), hjust = -0.15, size = 6 / .pt, colour = "grey15") +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000), labels = c("0.001", "0.01", "0.1", "1", "10", "100", "1,000"),
                limits = c(0.001, NA), expand = expansion(mult = c(0, 0.16))) +
  scale_y_discrete(limits = rev(d_lv)) + scale_colour_manual(values = d_col, guide = "none") +
  labs(x = "Projected CPU-hours on one core (log scale)", y = NULL, tag = "d") +
  theme_paper() + theme(panel.grid.major.y = element_blank(), plot.margin = margin(4, 14, 4, 4))
p <- pABC / free(pD) + plot_layout(heights = c(1.4, 1))

fmt <- function(v) format(v, big.mark = ",", scientific = FALSE, trim = TRUE)
reps <- unique(tm$count)
cap <- paste(
  "Figure S5. Computational cost of G-MoQ and direct mixture maximum likelihood under the two-component mixture model.",
  sprintf(paste("Panels a-c give the cost of one matched panel of %s simulated variants at three sample sizes, for the",
                "whole G-MoQ pipeline with the pooled or the 1,000-block jackknife covariance and for direct mixture",
                "maximum likelihood. **a**, Time to fit the estimates and their standard errors (s). **b**, Wall time of the",
                "whole process (s). **c**, Peak resident memory (MiB). Each G-MoQ time includes the one-off baseline mixture fit,",
                "Stage 1 and Stage 2 for the panel, so it is dominated by costs not repeated per variant (per-variant",
                "costs, Fig. 2c). Each point is the arithmetic mean of %s timing repeats of one method at one sample size,",
                "on log axes; no error bars are shown."),
          unique(tm$M), if (length(reps) == 1 && reps == 3) "three" else paste(reps, collapse = "/")),
  sprintf(paste("**d**, Projected cost of %s million variants at %s participants on one core, from the per-variant",
                "costs of Fig. 2c; bars start at 0.001 CPU-hours. Direct mixture maximum likelihood",
                "projects to %s CPU-hours; G-MoQ Stages 1 and 2, %.2f CPU-hours (%.0f times faster); Stage 2 alone given",
                "the Stage 1 output, %.1f seconds (about %s times faster; separately measured in-memory Stage 2 time 6.18 s). One-off per-trait costs are excluded."),
          G_snps / 1e6, fmt(4e5), fmt(round(g$cpu_h[1])), g$cpu_h[2], g$cpu_h[1] / g$cpu_h[2], g$cpu_h[3] * 3600,
          fmt(round(g$cpu_h[1] / g$cpu_h[3], -4))),
  "No statistical tests are shown.",
  "ML, maximum likelihood; SE, standard error; MiB, mebibyte; CPU, central processing unit.")
src <- rbind(tm[, .(panel = "abc", n, M, method, metric, value = mean_value, timing_repeats = count)],
             g[, .(panel = "d", n = 4e5, M = G_snps, method = as.character(what), metric = "genome-wide CPU-hours",
                   value = cpu_h, timing_repeats = NA_integer_)])
save_figure(p, "figS05_mixture_timing", src, cap, width = 183, height = 110,
            results_used = c("08_simulations/speed_benchmark", "08_simulations/speed_benchmark"))
