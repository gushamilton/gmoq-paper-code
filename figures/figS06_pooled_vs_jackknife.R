#!/usr/bin/env Rscript
# Supplementary Figure S6: pooled and jackknife covariance across sample size and allele frequency
#
# What: the pooled Stage 1 covariance assumes the genotype effect does not change the error
#   structure much; the 1,000-block jackknife makes no such assumption but is slower. On
#   simulated normal phenotypes (no effect, or a dispersion effect) at every combination of
#   sample size and minor allele frequency, the Stage 2 location and dispersion estimates are
#   checked against their projection targets under each covariance.
#   a: bias over the empirical s.d. (pooled and jackknife share point estimates).
#   b: mean SE over the empirical s.d., for each covariance.
#   c: coverage of nominal 95% confidence intervals, for each covariance.
# Paper: Supplementary Fig. S6
# Input: results/08_simulations/summary/table_s2_numerical_summary.csv
# Output: paper_assets/figS06_pooled_vs_jackknife.{png,pdf}, _source.csv, _caption.md
# Runs: locally from results/ in about 2 seconds

suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(data.table)})
source(file.path(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE), "figures", "theme.R"))

## 1. Load and label rows -----------------------------------------------------------------
sim <- as.data.table(read_result("08_simulations/summary", "table_s2_numerical_summary"))[study == "S6"]
sim[, cov_f := factor(method, levels = c("pooled", "jackknife"), labels = c("Pooled", "Jackknife"))]
cell_lab <- function(n, maf) sprintf("N %s, MAF %s", format(n, big.mark = ",", scientific = FALSE, trim = TRUE), maf)
sim[, cell := factor(cell_lab(n, maf), levels = rev(unique(cell_lab(n, maf)[order(n, -maf)])))]
# Facet rows: scenario x parameter, with Greek letters as plotmath.
group_lab <- c('paste("Null, ", beta[mu])', 'paste("Null, ", beta[sigma^2])',
               'paste("Dispersion, ", beta[mu])', 'paste("Dispersion, ", beta[sigma^2])')
sim[, grp := factor(group_lab[(scenario == "normal_dispersion") * 2 + (parameter == "beta_sigma2") + 1], levels = group_lab)]
sim[, mcse := empirical_sd / sqrt(success)]          # Monte Carlo standard error of the bias
cov_col <- c(Pooled = paper_palette[["gmoq"]], Jackknife = "#7fb3d5")
dodge <- position_dodge(width = 0.6)
row_theme <- theme_paper() + theme(panel.grid.major.y = element_blank(), strip.text.y.right = element_text(angle = 0, hjust = 0))

## 2. Panels a-c --------------------------------------------------------------------------
pA <- ggplot(sim[method == "pooled"], aes(bias / empirical_sd, cell)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_errorbar(aes(xmin = (bias - 1.96 * mcse) / empirical_sd, xmax = (bias + 1.96 * mcse) / empirical_sd), width = 0,
                linewidth = 0.55, colour = "grey30", orientation = "y") +
  geom_point(size = 1.5, colour = "grey20") +
  facet_grid(grp ~ ., labeller = label_parsed) + labs(x = "Bias / empirical SD", y = NULL, tag = "a") +
  row_theme + theme(strip.text.y.right = element_blank())
pB <- ggplot(sim, aes(mean_se / empirical_sd, cell, colour = cov_f)) +
  geom_vline(xintercept = 1, linewidth = 0.3, linetype = 2, colour = "grey45") +
  geom_point(size = 1.5, position = dodge) +
  facet_grid(grp ~ ., labeller = label_parsed) + scale_colour_manual(values = cov_col, name = NULL) +
  labs(x = "Mean SE / empirical SD", y = NULL, tag = "b") +
  row_theme + theme(axis.text.y = element_blank(), strip.text.y.right = element_blank(), legend.position = "top")
pC <- ggplot(sim, aes(coverage, cell, colour = cov_f)) +
  geom_vline(xintercept = 0.95, linewidth = 0.3, linetype = 2, colour = "grey45") +
  geom_errorbar(aes(xmin = coverage - 1.96 * coverage_mcse, xmax = coverage + 1.96 * coverage_mcse), width = 0,
                linewidth = 0.5, position = dodge, orientation = "y") +
  geom_point(size = 1.5, position = dodge) +
  facet_grid(grp ~ ., labeller = label_parsed) + scale_colour_manual(values = cov_col, guide = "none") +
  labs(x = "95% CI coverage", y = NULL, tag = "c") +
  row_theme + theme(axis.text.y = element_blank())

## 3. Assemble and save ---------------------------------------------------------------------
fig <- (pA | pB | pC) + plot_layout(widths = c(1.1, 1, 1.15))
fmt <- function(v) format(v, big.mark = ",", scientific = FALSE, trim = TRUE)
or3 <- function(v) paste(paste(head(v, -1), collapse = ", "), "or", tail(v, 1))
cap <- paste(
  "Figure S6. Pooled and 1,000-block jackknife covariance across sample size and minor allele frequency.",
  sprintf(paste("Simulated normal phenotypes with no genetic effect (Null) or a genetic effect on dispersion (Dispersion),",
                "in every combination (rows) of %s participants and minor allele frequency %s; Stage 2 location (β_μ) and dispersion",
                "(β_σ²) estimates by ordinary least squares, against their population projection targets (Methods)."),
          or3(fmt(sort(unique(sim$n)))), or3(sort(unique(sim$maf), decreasing = TRUE))),
  sprintf("Each point summarizes %s replicates of one row; MC interval, ±1.96 Monte Carlo standard errors.", fmt(unique(sim$success))),
  "**a**, Bias as a fraction of the empirical standard deviation of the estimates, with 95% MC intervals;",
  "pooled and jackknife share point estimates.",
  "**b**, Mean SE over the empirical standard deviation, for each covariance; the dashed line marks 1.",
  "**c**, Coverage of the nominal 95% confidence interval for each covariance, with 95% MC intervals;",
  "the dashed line marks 0.95.",
  "No statistical tests are shown.",
  "MC, Monte Carlo; SD, standard deviation; SE, standard error; CI, confidence interval; MAF, minor allele frequency.")
save_figure(fig, "figS06_pooled_vs_jackknife", sim[, .(scenario, parameter, n, maf, method, bias, empirical_sd, mean_se, coverage,
            coverage_mcse, replicates = success)], cap, width = 183, height = 190,
            results_used = "08_simulations/summary")
