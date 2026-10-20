#!/usr/bin/env Rscript
# Supplementary Figure S14: interaction of BMI leads with the BMI polygenic score, by scale of BMI
#
# What: the common relative interaction kappa = delta / beta (interaction with the polygenic
#   score per s.d., over main effect) across the 639 BMI leads depends on the scale BMI is
#   measured on. It is shown on raw BMI, log BMI, the rank-inverse-normal scale and the implied
#   natural scale (under which the leads' common profile is a pure shift).
#   Single panel: kappa per scale with 95% bootstrap confidence intervals over leads.
# Paper: Supplementary Fig. S14
# Input: results/05_bmi/prs_interaction/kappa_by_scale.csv
# Output: paper_assets/figS14_prs_interaction_by_scale.{png,pdf}, _source.csv, _caption.md
# Runs: locally from results/ in about 2 seconds

suppressPackageStartupMessages({library(ggplot2); library(data.table)})
source(file.path(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE), "figures", "theme.R"))

## 1. Load (per-lead regressions with heteroskedasticity-robust HC0 standard errors) --------
kappa <- as.data.table(read_result("05_bmi/prs_interaction", "kappa_by_scale"))[se_type == "HC0"]
scale_lab <- c(raw = "Raw BMI", log = "Log BMI", int = "Rank-inverse\nnormal", natural = "Implied\nnatural")
stopifnot(setequal(kappa$scale, names(scale_lab)), uniqueN(kappa$n_leads) == 1)
kappa[, scale_f := factor(scale_lab[scale], levels = scale_lab)]

## 2. Plot ----------------------------------------------------------------------------------
fig <- ggplot(kappa, aes(scale_f, kappa)) +
  geom_col(fill = paper_palette[["shift"]], width = 0.6, alpha = 0.85) +
  geom_errorbar(aes(ymin = kappa_lo, ymax = kappa_hi), width = 0.2, linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.3f", kappa), y = kappa_hi), vjust = -0.5, size = 6 / .pt) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = expression(paste("Relative interaction ", kappa, " (per s.d. of PRS)"))) +
  theme_paper()

## 3. Save ----------------------------------------------------------------------------------
src <- kappa[, .(scale, n_leads, kappa, kappa_se_boot, kappa_lo, kappa_hi)]
k <- function(s) kappa[scale == s]
cap <- paste(
  "Figure S14. Interaction of BMI lead variants with the BMI polygenic score on four scales of BMI.",
  sprintf("Each bar is the common relative interaction κ = δ/β across the %d BMI lead variants (BMI-increasing allele;",
          kappa$n_leads[1]),
  "n = 461,320), where β is a variant's main effect and δ its interaction with a polygenic risk score (PRS) per standard deviation",
  "of the score, estimated with a score that excludes the variant's chromosome, on raw BMI, log BMI, rank-inverse-normal",
  "BMI and the implied natural scale, under which the leads' common quantile-effect profile is a pure shift (Methods).",
  "Error bars are 95% bootstrap confidence intervals over leads; per-lead regressions used heteroskedasticity-robust",
  sprintf("standard errors. κ is %.3f on raw BMI, %.3f on log BMI, %.3f on the rank-inverse-normal scale and %.3f on the",
          k("raw")$kappa, k("log")$kappa, k("int")$kappa, k("natural")$kappa),
  "implied natural scale. No statistical tests are shown.",
  "BMI, body mass index; PRS, polygenic risk score; s.d., standard deviation.")
save_figure(fig, "figS14_prs_interaction_by_scale", src, cap, width = 89, height = 70,
            results_used = "05_bmi/prs_interaction")
