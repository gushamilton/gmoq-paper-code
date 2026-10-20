#!/usr/bin/env Rscript
# Supplementary Figure S3: Stage 2 estimates against their population projection targets
#
# What: when the model does not contain the generating effect exactly, Stage 2 estimates the
#   projection of that effect onto the model's weight columns. In simulated data (normal
#   phenotypes with location and/or dispersion effects, and a two-component mixture), the
#   estimates are compared with that target, with pooled and jackknife covariance.
#   a: bias over RMSE, with 95% Monte Carlo intervals.
#   b: mean SE over RMSE.
#   c: coverage of nominal 95% confidence intervals.
# Paper: Supplementary Fig. S3
# Input: results/08_simulations/summary/table_s2_numerical_summary.csv
# Output: paper_assets/figS03_projection_targets.{png,pdf}, _source.csv, _caption.md
# Runs: locally from results/ in about 2 seconds

suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(data.table)})
source(file.path(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE), "figures", "theme.R"))

## 1. Load and label rows -----------------------------------------------------------------
sim <- as.data.table(read_result("08_simulations/summary", "table_s2_numerical_summary"))[study == "S3"]
scen_lab <- c(normal_location = "Normal, location", normal_dispersion = "Normal, dispersion", normal_both = "Normal, both",
              mixture_left = "Mixture, left")
# Greek letters as plotmath, which every PDF device can draw.
par_lab <- c(beta_mu = "beta[mu]", beta_sigma2 = "beta[sigma^2]", beta_1 = "beta[1]")
sim[, row_s := sprintf('paste("%s: ", %s)', scen_lab[scenario], par_lab[parameter])]
sim[, row := factor(row_s, levels = rev(unique(row_s[order(match(scenario, names(scen_lab)), parameter)])))]
sim[, cov_f := factor(covariance_source, levels = c("pooled", "jackknife"), labels = c("Pooled", "Jackknife"))]
sim[, mcse_bias := rmse / sqrt(success)]          # Monte Carlo standard error of the bias
cov_col <- c(Pooled = paper_palette[["gmoq"]], Jackknife = "#7fb3d5")
dodge <- position_dodge(width = 0.55)
common <- function(p) p + scale_colour_manual(values = cov_col, name = NULL) +
  scale_y_discrete(labels = function(b) parse(text = b)) + labs(y = NULL) +
  theme_paper() + theme(panel.grid.major.y = element_blank())

## 2. Panels a-c --------------------------------------------------------------------------
pA <- common(ggplot(sim, aes(bias / rmse, row, colour = cov_f)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_errorbar(aes(xmin = (bias - 1.96 * mcse_bias) / rmse, xmax = (bias + 1.96 * mcse_bias) / rmse), width = 0,
                linewidth = 0.6, position = dodge, orientation = "y") +
  geom_point(size = 1.6, position = dodge) +
  labs(x = "Bias / RMSE", tag = "a")) + theme(legend.position = "top", legend.justification = "left")
pB <- common(ggplot(sim, aes(mean_se / rmse, row, colour = cov_f)) +
  geom_vline(xintercept = 1, linewidth = 0.3, linetype = 2, colour = "grey45") +
  geom_point(size = 1.6, position = dodge) +
  labs(x = "Mean SE / RMSE", tag = "b")) + theme(axis.text.y = element_blank(), legend.position = "none")
pC <- common(ggplot(sim, aes(coverage, row, colour = cov_f)) +
  geom_vline(xintercept = 0.95, linewidth = 0.3, linetype = 2, colour = "grey45") +
  geom_errorbar(aes(xmin = coverage - 1.96 * coverage_mcse, xmax = coverage + 1.96 * coverage_mcse), width = 0,
                linewidth = 0.6, position = dodge, orientation = "y") +
  geom_point(size = 1.6, position = dodge) +
  labs(x = "95% CI coverage", tag = "c")) + theme(axis.text.y = element_blank(), legend.position = "none")

## 3. Assemble and save ---------------------------------------------------------------------
fig <- (pA | pB | pC) + plot_layout(widths = c(1.2, 1, 1))
fmt <- function(v) format(v, big.mark = ",", scientific = FALSE, trim = TRUE)
cap <- paste(
  "Figure S3. Accuracy and calibration of Stage 2 estimates against their population projection targets.",
  sprintf(paste("Simulated data on %s participants under five scenario-parameter combinations (rows): a normal",
                "phenotype with a genetic effect on location, on dispersion or on both, with location (β_μ) and",
                "dispersion (β_σ²) estimates, and a two-component mixture with an effect on the left component (β_1)."),
          fmt(unique(sim$n))),
  "The target is the projection of the generating effect onto the model's weight columns (Methods).",
  sprintf(paste("Each point summarizes %s replicates of one row with the pooled or the 1,000-block jackknife covariance;",
                "the two share point estimates. MC interval, ±1.96 Monte Carlo standard errors."), fmt(unique(sim$success))),
  "**a**, Bias as a fraction of the RMSE, with 95% MC intervals. **b**, Mean SE over the RMSE; the dashed line marks 1.",
  "**c**, Coverage of the nominal 95% confidence interval, with 95% MC intervals; the dashed line marks 0.95.",
  "No statistical tests are shown.",
  "MC, Monte Carlo; RMSE, root mean square error; SE, standard error; CI, confidence interval.")
save_figure(fig, "figS03_projection_targets", sim[, .(scenario, parameter, covariance_source, n, bias, rmse, mean_se, coverage,
            coverage_mcse, replicates = success)], cap, width = 183, height = 65,
            results_used = "08_simulations/summary")
