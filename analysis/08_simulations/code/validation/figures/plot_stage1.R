#!/usr/bin/env Rscript
# Draft plots of the Stage 1 comparison
#
# What: draws working plots of the coefficient agreement, standard errors and timings exported
#   by export_stage1.py. These are run-time checks only; the paper's Supplementary Figs. S1 and
#   S2 are drawn by figures/figS01_stage1_identity.R and figS02_stage1_timing.R from the
#   imported summaries.
# Paper: none directly (studies 1 and 2)
# Input: the tables directory written by export_stage1.py
# Output: PNG/PDF working plots in a new figure directory
# Runs: on a cluster, called at the end of run_stage1.py, under a minute

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript plot_stage1.R TABLE_DIR FIGURE_DIR", call. = FALSE)
}
table_dir <- normalizePath(args[[1L]], mustWork = TRUE)
figure_dir <- normalizePath(args[[2L]], mustWork = FALSE)
if (file.exists(figure_dir)) {
  stop("refusing existing exclusive figure directory: ", figure_dir, call. = FALSE)
}
host <- unname(Sys.info()[["nodename"]])
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("ggplot2 is required in the R 4.5.1 module", call. = FALSE)
}
suppressPackageStartupMessages(library(ggplot2))
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

need <- function(x, columns, label) {
  if (!all(columns %in% names(x))) stop(label, " is missing required columns", call. = FALSE)
  x
}
finite_columns <- function(x, columns, label) {
  for (column in columns) {
    value <- suppressWarnings(as.numeric(x[[column]]))
    if (any(!is.finite(value))) stop(label, " has non-finite ", column, call. = FALSE)
    x[[column]] <- value
  }
  x
}
s1 <- read.csv(file.path(table_dir, "S1_pairs.csv"), stringsAsFactors = FALSE,
               check.names = FALSE)
s1 <- need(s1, c("snp_id", "tau", "beta_gmoq", "beta_rif_ols", "difference"), "S1_pairs.csv")
s1 <- finite_columns(s1, c("tau", "beta_gmoq", "beta_rif_ols", "difference"), "S1_pairs.csv")
s2 <- read.csv(file.path(table_dir, "S2_timings.csv"), stringsAsFactors = FALSE,
               check.names = FALSE)
s2 <- need(s2, c("workload", "method", "repeat", "n", "m", "fit_seconds",
                 "wall_seconds", "peak_rss_kib"), "S2_timings.csv")
s2 <- finite_columns(s2, c("repeat", "n", "m", "fit_seconds", "wall_seconds",
                            "peak_rss_kib"), "S2_timings.csv")
if (nrow(s2) != 36L || !all(s2[["repeat"]] %in% 1:3)) stop("S2 timings must contain 36 three-repeat rows", call. = FALSE)

save_plot <- function(plot, stem, width, height) {
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot, device = "pdf",
         width = width, height = height, units = "in")
  ggsave(file.path(figure_dir, paste0(stem, ".png")), plot, device = "png",
         width = width, height = height, units = "in", dpi = 300)
}

p_s1 <- ggplot(s1, aes(x = beta_rif_ols, y = beta_gmoq, colour = tau)) +
  geom_point(size = 1.25, alpha = 0.72) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.45, linetype = 2) +
  scale_colour_gradient(low = "#2166ac", high = "#b2182b") +
  labs(title = "Numerical equivalence to ordinary RIF-OLS",
       subtitle = "Largest matched panel; identity line shown",
       x = "Ordinary RIF-OLS coefficient (reference)",
       y = "G-MoQ Stage 1 coefficient", colour = "tau") +
  theme_bw(base_size = 11) + theme(legend.position = "right")
save_plot(p_s1, "S1", 7.2, 5.5)

p_diff <- ggplot(s1, aes(x = tau, y = difference, colour = tau)) +
  geom_hline(yintercept = 0, linewidth = 0.45, linetype = 2) +
  geom_point(size = 1.15, alpha = 0.72) +
  scale_colour_gradient(low = "#2166ac", high = "#b2182b") +
  labs(title = "Numerical difference from ordinary RIF-OLS",
       subtitle = "G-MoQ coefficient minus ordinary RIF-OLS coefficient",
       x = "tau", y = "Difference", colour = "tau") +
  theme_bw(base_size = 11)
save_plot(p_diff, "S1_difference", 7.2, 5.5)

s2$method_label <- unname(c(gmoq="G-MoQ pooled", ordinary_rif_ols="ordinary RIF-OLS", gmoq_jackknife="G-MoQ jackknife")[s2$method])
workload_info <- unique(s2[c("workload", "n", "m")])
workload_info <- workload_info[order(workload_info$n, workload_info$m, workload_info$workload), ]
s2$workload <- factor(s2$workload, levels = workload_info$workload)
workload_labels <- setNames(
  paste0("N=", formatC(as.integer(workload_info$n), format="d", big.mark=","),
         "\nM=", formatC(as.integer(workload_info$m), format="d", big.mark=",")),
  workload_info$workload)
long <- rbind(
  data.frame(workload = s2$workload, method = s2$method_label, replicate_id = s2[["repeat"]],
             metric = "fit-only seconds", value = s2$fit_seconds),
  data.frame(workload = s2$workload, method = s2$method_label, replicate_id = s2[["repeat"]],
             metric = "whole-process seconds", value = s2$wall_seconds),
  data.frame(workload = s2$workload, method = s2$method_label, replicate_id = s2[["repeat"]],
             metric = "peak RSS MiB", value = s2$peak_rss_kib / 1024)
)
long$metric <- factor(long$metric,
  levels = c("fit-only seconds", "whole-process seconds", "peak RSS MiB"))
long$method <- factor(long$method, levels = c("G-MoQ pooled", "ordinary RIF-OLS", "G-MoQ jackknife"))
p_s2 <- ggplot(long[long$method != "G-MoQ jackknife", ], aes(x = workload, y = value, colour = method, group = method,
                          shape = factor(replicate_id))) +
  geom_point(position = position_jitterdodge(jitter.width = 0.06, dodge.width = 0.68, seed = 940200),
             size = 2.2, alpha = 0.88) +
  stat_summary(aes(group = method), fun = median, geom = "point", shape = 95,
               size = 8, position = position_dodge(width = 0.68), show.legend = FALSE) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  scale_x_discrete(labels = workload_labels) +
  scale_colour_manual(values = c("G-MoQ pooled" = "#2166ac", "ordinary RIF-OLS" = "#b2182b", "G-MoQ jackknife" = "#008837")) +
  labs(title = "Stage 1 coefficient and SE timing and memory measurements",
       subtitle = "All three repeats shown; thick marks are medians; common RIF setup is included in whole-process time",
       x = "Matched workload (isolated M=1 versus M-panel scan)", y = NULL,
       colour = "Method", shape = "Repeat") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(size = 8), legend.position = "bottom")
save_plot(p_s2, "S2", 10.0, 4.8)
p_jk <- (p_s2 %+% long[long$method != "G-MoQ pooled", ]) +
  labs(title = "Stage 1 coefficient and SE timing: jackknife sensitivity")
save_plot(p_jk, "S2_jackknife", 10.0, 4.8)
se <- read.csv(file.path(table_dir, "S1_SE.csv"))
se_long <- rbind(data.frame(reference=se$se_ordinary_rif_ols, gmoq=se$se_gmoq, source="Pooled"),
                 data.frame(reference=se$se_ordinary_rif_ols, gmoq=se$se_gmoq_jackknife, source="Jackknife"))
p_se <- ggplot(se_long, aes(reference, gmoq)) + geom_point(size=.5, alpha=.3) +
  geom_abline(slope=1, intercept=0, linetype=2) + facet_wrap(~source) + theme_bw() +
  labs(title="Stage 1 standard errors", subtitle="Covariance conventions differ; identity line is a reference, not an equivalence requirement",
       x="Conventional RIF-OLS SE", y="G-MoQ SE")
save_plot(p_se, "S1_SE", 9, 4.5)

metadata <- c(
  "Stage 1 figure rendering metadata", paste0("table_dir=", table_dir),
  paste0("figure_dir=", figure_dir), paste0("R=", R.version.string),
  paste0("ggplot2=", as.character(packageVersion("ggplot2"))),
  paste0("script=", normalizePath(commandArgs()[1L], mustWork = FALSE)),
  "Plots use ggplot2 and base R only.",
  "Timings include coefficients and SEs; pooled, conventional OLS and jackknife covariance conventions differ.",
  "The ordinary baseline uses np.linalg.lstsq(design, RIF45responses) once per SNP."
)
writeLines(metadata, file.path(figure_dir, "plot_metadata.txt"))
