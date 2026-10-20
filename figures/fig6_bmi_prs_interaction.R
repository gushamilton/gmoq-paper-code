#!/usr/bin/env Rscript
# Figure 6: BMI variant effects grow with polygenic background
#
# What: the common tilt of the BMI lead profiles (Fig. 4c) means that alleles have larger
#   effects in people with higher BMI. One explanation is that each allele's effect grows
#   with the rest of the genome. For each of the 639 leads, log BMI was regressed on the
#   variant, a polygenic score (PRS) that excludes the variant's chromosome, and their product.
#   a: distribution of the per-lead interaction z against the standard normal null.
#   b: interaction against main effect per lead, with the common relative interaction kappa.
#   c: effect of the FTO lead rs1421085 by PRS decile.
#   d: the observed common tilt beside the tilt predicted by four candidate modifiers.
# Paper: Fig. 6
# Input: results/05_bmi/prs_interaction/{per_lead_prs,kappa_by_scale,fto_prs_deciles,
#   fig5d_curves,fig5d_summary}.csv
# Output: paper_assets/fig6_bmi_prs_interaction.{png,pdf}, _source.csv, _caption.md
# Runs: locally from results/ in about 2 seconds

suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(data.table)})
source(file.path(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE), "figures", "theme.R"))

## 1. Load ----------------------------------------------------------------------------------
# Per-lead regressions on log BMI, BMI-increasing allele, heteroscedasticity-robust (HC0) SEs.
leads <- as.data.table(read_result("05_bmi/prs_interaction", "per_lead_prs"))
kappa_log <- as.data.table(read_result("05_bmi/prs_interaction", "kappa_by_scale"))[scale == "log" & se_type == "HC0"]
fto <- as.data.table(read_result("05_bmi/prs_interaction", "fto_prs_deciles"))[scale == "log"]
curves <- as.data.table(read_result("05_bmi/prs_interaction", "fig5d_curves"))
curve_fit <- as.data.table(read_result("05_bmi/prs_interaction", "fig5d_summary"))
fto_contrast <- as.data.table(read_result("05_bmi/prs_interaction", "fto_prs_contrast"))[scale == "log"]

## 2. Panel a: interaction z against the null ---------------------------------------------
leads[, z_int := gamma_log / se_gamma_log_hc0]
null_curve <- data.table(z = seq(-4, 4, length.out = 300)); null_curve[, d := dnorm(z)]
pA <- ggplot(leads, aes(z_int)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.25, boundary = 0, fill = "#5f9ea0", colour = "white",
                 linewidth = 0.15) +
  geom_line(data = null_curve, aes(z, d), colour = "grey15", linewidth = 0.7) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3, colour = "grey45") +
  labs(x = "SNP x PRS interaction z (BMI-increasing allele, HC0)", y = "Density", tag = "a") +
  theme_paper()

## 3. Panel b: main effect against interaction --------------------------------------------
pB <- ggplot(leads, aes(beta_log, gamma_log)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey55") +
  geom_errorbar(aes(ymin = gamma_log - 1.96 * se_gamma_log_hc0, ymax = gamma_log + 1.96 * se_gamma_log_hc0),
                width = 0, colour = "grey85", linewidth = 0.2) +
  geom_point(size = 0.7, alpha = 0.7, colour = paper_palette[["shift"]], stroke = 0) +
  geom_abline(slope = kappa_log$kappa, intercept = 0, colour = paper_palette[["spread"]], linewidth = 0.7) +
  labs(x = "SNP main effect beta_i (log BMI)", y = "SNP x PRS interaction gamma_i (per PRS SD)", tag = "b") +
  theme_paper()

## 4. Panel c: FTO by PRS decile -------------------------------------------------------------
fto[, extreme := prs_decile %in% c(1, 10)]
pC <- ggplot(fto, aes(prs_decile, beta)) +
  geom_errorbar(aes(ymin = beta - 1.96 * se_hc0, ymax = beta + 1.96 * se_hc0, colour = extreme), width = 0, linewidth = 0.6) +
  geom_point(aes(colour = extreme), size = 1.9) +
  scale_colour_manual(values = c(`TRUE` = paper_palette[["highlight"]], `FALSE` = "grey35"), guide = "none") +
  scale_x_continuous(breaks = 1:10) +
  labs(x = "PRS decile (chromosome 16 excluded)", y = "FTO rs1421085 effect on log BMI (95% CI)", tag = "c") +
  theme_paper()

## 5. Panel d: observed and predicted common profiles -------------------------------------
# b(tau)/beta is the common lead profile per unit mean effect; each prediction assumes the
# allele's effect grows with a modifier correlated with BMI, and F is its fit to the observed.
curve_lab <- c(observed = "observed (95% band)",
               pred_bmi = sprintf("BMI itself (F = %.2f)", curve_fit[hyp_key == "bmi", F]),
               pred_twin = sprintf("twin h2 0.5 (F = %.2f)", curve_fit[hyp_key == "twin", F]),
               pred_h2snp = sprintf("h2_SNP %.3f (F = %.2f)", curve_fit[hyp_key == "h2snp", omega], curve_fit[hyp_key == "h2snp", F]),
               pred_prs = sprintf("measured PRS (F = %.2f)", curve_fit[hyp_key == "prs", F]))
curves_long <- melt(curves, id.vars = "tau", measure.vars = names(curve_lab), variable.name = "series", value.name = "value")
curves_s <- smooth_curve(curves_long, "tau", "value", by = "series")
curves_s[, series := factor(series, levels = names(curve_lab))]
band <- merge(smooth_curve(curves, "tau", "observed_lo")[, .(x, lo = y)], smooth_curve(curves, "tau", "observed_hi")[, .(x, hi = y)], by = "x")
curve_col <- c(observed = "black", pred_bmi = paper_palette[["spread"]], pred_twin = paper_palette[["bend"]],
               pred_h2snp = paper_palette[["shift"]], pred_prs = paper_palette[["mixture"]])
pD <- ggplot(curves_s, aes(x, y, colour = series)) +
  geom_ribbon(data = band, aes(x, ymin = lo, ymax = hi), inherit.aes = FALSE, fill = "grey75", alpha = 0.5) +
  geom_hline(yintercept = 1, linetype = 2, linewidth = 0.3, colour = "grey50") +
  geom_line(aes(linewidth = series == "observed"), lineend = "round") +
  scale_linewidth_manual(values = c(`TRUE` = 0.9, `FALSE` = 0.65), guide = "none") +
  scale_colour_manual(values = curve_col, labels = curve_lab, name = NULL) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", ".25", ".5", ".75", "1")) +
  labs(x = "Quantile level tau", y = "b(tau) / beta", tag = "d") +
  theme_paper() + legend_inside(0.02, 0.98, c(0, 1))

## 6. Assemble and save ---------------------------------------------------------------------
fig <- (pA | pB) / (pC | pD) + plot_layout(heights = c(1, 1))
src <- rbind(
  leads[, .(panel = "A/B", series = snp_id, x = beta_log, y = gamma_log, se = se_gamma_log_hc0, z = z_int)],
  fto[, .(panel = "C", series = paste0("decile_", prs_decile), x = prs_decile, y = beta, se = se_hc0, z = NA_real_)],
  curves_long[, .(panel = "D", series = as.character(series), x = tau, y = value, se = NA_real_, z = NA_real_)])
cap <- paste(
  "Figure 6. BMI variant effects grow with polygenic background.",
  "Log BMI, 639 lead variants aligned to the BMI-increasing allele; n = 461,320. For each lead, log BMI was regressed",
  "on the variant, a polygenic score excluding that variant's chromosome, their interaction, age, sex and ten genetic",
  "principal components, with heteroscedasticity-robust standard errors. **a**, Distribution of the interaction z",
  sprintf("statistic against the standard normal null; interactions are positive at %.0f%% of leads. **b**, Interaction against",
          100 * mean(leads$gamma_log > 0)),
  "main effect per lead, with 95% confidence intervals; the red line is the common relative interaction",
  sprintf("kappa = %.3f (95%% confidence interval %.3f to %.3f). **c**, Effect of the FTO lead rs1421085 on log BMI by",
          kappa_log$kappa, kappa_log$kappa_lo, kappa_log$kappa_hi),
  sprintf("decile of the polygenic score, with 95%% confidence intervals; the top decile effect is %.2f times the bottom",
          fto_contrast$ratio_top_over_bottom),
  sprintf("(95%% confidence interval %.2f to %.2f). **d**,", fto_contrast$ratio_lo, fto_contrast$ratio_hi),
  "The common tilt of the lead profiles, observed and as predicted if",
  "each allele's effect grew with a modifier correlated with BMI, for four candidate modifiers. Amplification by",
  "polygenic background is one reading of these results; a convex measurement scale and environments correlated with",
  "genetic liability predict the same pattern.")
save_figure(fig, "fig6_bmi_prs_interaction", src, cap, width = 183, height = 150,
            results_used = "05_bmi/prs_interaction")
