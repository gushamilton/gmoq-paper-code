#!/usr/bin/env Rscript
# Supplementary Figure S1: Stage 1 reproduces independently fitted RIF-OLS regressions
#
# What: Stage 1 fits all 45 quantile levels at once by RIF-OLS (ordinary least squares on the
#   recentred influence function). In simulated data its coefficients and standard errors
#   are compared with separately fitted ordinary RIF-OLS regressions, coefficient by
#   coefficient. Only the aggregate summary is plotted.
#   a: coefficient difference, mean and +-2 s.d. across coefficients.
#   b: ratio of the Stage 1 SE (pooled and jackknife covariance) to the ordinary SE.
# Paper: Supplementary Fig. S1
# Input: results/08_simulations/summary/{table_s2_numerical_summary,table_s1_design_provenance}.csv
# Output: paper_assets/figS01_stage1_identity.{png,pdf}, _source.csv, _caption.md
# Runs: locally from results/ in about 2 seconds

suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(data.table)})
source(file.path(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE), "figures", "theme.R"))

## 1. Load ----------------------------------------------------------------------------------
sim <- as.data.table(read_result("08_simulations/summary", "table_s2_numerical_summary"))[study == "S1"]
design <- as.data.table(read_result("08_simulations/summary", "table_s1_design_provenance"))[study == "S1"]
coef_diff <- sim[quantity == "coefficient_minus_RIF_OLS"]
se_ratio <- sim[quantity == "SE_ratio_to_ordinary"]
se_ratio[, method_f := factor(method, levels = c("pooled", "jackknife"), labels = c("Pooled", "Jackknife"))]
se_ratio[, `:=`(mean = 1 + bias, sd = empirical_sd)]

## 2. Panel a: coefficient difference (in units of 1e-14) --------------------------------
pA <- ggplot(coef_diff, aes(x = 1)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_errorbar(aes(ymin = (bias - 2 * empirical_sd) * 1e14, ymax = (bias + 2 * empirical_sd) * 1e14), width = 0,
                linewidth = 0.7, colour = paper_palette[["gmoq"]]) +
  geom_point(aes(y = bias * 1e14), size = 2.2, colour = paper_palette[["gmoq"]]) +
  scale_x_continuous(limits = c(0.7, 1.3), breaks = NULL) +
  labs(x = NULL, y = expression("Coefficient difference (" * 10^-14 * ")"), tag = "a") +
  theme_paper()

## 3. Panel b: SE ratio -------------------------------------------------------------------
pB <- ggplot(se_ratio, aes(mean, method_f, colour = method_f)) +
  geom_vline(xintercept = 1, linewidth = 0.3, linetype = 2, colour = "grey45") +
  geom_errorbar(aes(xmin = mean - 2 * sd, xmax = mean + 2 * sd), width = 0, linewidth = 0.8, orientation = "y") +
  geom_point(size = 2.2) +
  scale_colour_manual(values = c(Pooled = paper_palette[["gmoq"]], Jackknife = "#7fb3d5"), guide = "none") +
  scale_x_continuous(limits = c(0.94, 1.06)) +
  labs(x = "SE ratio, G-MoQ / ordinary RIF-OLS", y = NULL, tag = "b") +
  theme_paper() + theme(panel.grid.major.y = element_blank())

## 4. Assemble and save ---------------------------------------------------------------------
fig <- (pA | pB) + plot_layout(widths = c(0.5, 2))
src <- sim[, .(quantity, method, bias, rmse, empirical_sd, se_agreement, attempted, success)]
# a x 10^b with Unicode superscripts, for the legend text
sci <- function(v, d = 1) {
  e <- floor(log10(abs(v))); m <- v / 10^e
  sprintf("%.*f × 10%s", d, m, chartr("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹", as.character(e)))
}
n_coef <- coef_diff$success; n_tau <- 45L
cap <- paste(
  "Figure S1. Identity of Stage 1 G-MoQ and independently fitted ordinary RIF-OLS regressions.",
  sprintf("Simulated data on %s participants; %s coefficients (%s variants × %d quantile levels), each fitted",
          format(as.numeric(design$N), big.mark = ",", scientific = FALSE), format(n_coef, big.mark = ","),
          format(n_coef / n_tau, big.mark = ","), n_tau),
  "by G-MoQ Stage 1 and by a separate ordinary RIF-OLS regression (Methods).",
  sprintf(paste("**a**, Difference between the G-MoQ and the ordinary RIF-OLS coefficient: point, mean over the %s",
                "coefficients (%s); bar, ±2 s.d. across coefficients (root mean square difference %s)."),
          format(n_coef, big.mark = ","), sci(coef_diff$bias), sci(coef_diff$rmse)),
  sprintf(paste("**b**, Ratio of the G-MoQ standard error to the ordinary RIF-OLS standard error, with the pooled",
                "covariance and with the 1,000-block jackknife: point, mean over the same %s coefficients; bar, ±2 s.d.",
                "Pooled, mean %.6f (s.d. %s; bar hidden by its point); jackknife, mean %.6f (s.d. %.4f). The dashed line marks a ratio of 1."),
          format(n_coef, big.mark = ","), se_ratio[method == "pooled", mean], sci(se_ratio[method == "pooled", sd]),
          se_ratio[method == "jackknife", mean], se_ratio[method == "jackknife", sd]),
  "No statistical tests are shown.",
  "RIF-OLS, recentred influence function ordinary least squares; SE, standard error; s.d., standard deviation.")
save_figure(fig, "figS01_stage1_identity", src, cap, width = 120, height = 55,
            results_used = "08_simulations/summary")
