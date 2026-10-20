#!/usr/bin/env Rscript
# Supplementary Figure S4: G-MoQ Stage 2 against direct mixture maximum likelihood
#
# What: on simulated two-component mixtures, the G-MoQ mixture estimates (Stage 1 then Stage 2,
#   pooled or jackknife covariance) are compared with refitting the mixture by maximum
#   likelihood at every variant, for the lower and upper component locations (delta_1,
#   delta_2) and membership (gamma).
#   a: bias by sample size, with 95% Monte Carlo intervals.
#   b: mean SE over RMSE.
#   c: coverage of 95% intervals in a second simulation, by mixture design (separation).
#   d: median time to fit one 100-variant panel against sample size.
# Paper: Supplementary Fig. S4
# Input: results/08_simulations/summary/{table_s2_numerical_summary,lowsep_sim_summary,
#   mixture_timings}.csv
# Output: paper_assets/figS04_mixture_recovery.{png,pdf}, _source.csv, _caption.md
# Runs: locally from results/ in about 2 seconds

suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(data.table)})
source(file.path(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE), "figures", "theme.R"))

method_lab <- c(gmoq_pooled = "G-MoQ, pooled", gmoq_jackknife = "G-MoQ, jackknife", mixture_ml = "Mixture ML")
method_col <- setNames(c(paper_palette[["gmoq"]], "#7fb3d5", paper_palette[["mixture"]]), method_lab)
# Greek letters as plotmath, which every PDF device can draw.
par_lab <- c(delta1 = "delta[1]", delta2 = "delta[2]", gamma = "gamma")
fmt <- function(v) format(v, big.mark = ",", scientific = FALSE, trim = TRUE)

## 1. Panels a, b: bias and SE calibration ------------------------------------------------
sim <- as.data.table(read_result("08_simulations/summary", "table_s2_numerical_summary"))[study == "S4"]
sim[, method_f := factor(method_lab[method], levels = method_lab)]
sim[, n_f := factor(sprintf("N = %s", fmt(n)), levels = sprintf("N = %s", fmt(sort(unique(n)))))]
sim[, par_f := factor(par_lab[parameter], levels = par_lab)]
sim[, mcse := rmse / sqrt(success)]               # Monte Carlo standard error of the bias
dodge <- position_dodge(width = 0.6)
pA <- ggplot(sim, aes(bias, n_f, colour = method_f)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_errorbar(aes(xmin = bias - 1.96 * mcse, xmax = bias + 1.96 * mcse), width = 0, linewidth = 0.6, position = dodge,
                orientation = "y") +
  geom_point(size = 1.7, position = dodge) +
  facet_wrap(~par_f, ncol = 1, strip.position = "left", labeller = label_parsed) +
  scale_colour_manual(values = method_col, name = NULL) +
  labs(x = "Bias", y = NULL, tag = "a") +
  theme_paper() + theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1, size = rel(1.15)),
                        panel.grid.major.y = element_blank(), legend.position = "none")
pB <- ggplot(sim, aes(mean_se / rmse, n_f, colour = method_f)) +
  geom_vline(xintercept = 1, linewidth = 0.3, linetype = 2, colour = "grey45") +
  geom_point(size = 1.7, position = dodge) +
  facet_wrap(~par_f, ncol = 1, strip.position = "left", labeller = label_parsed) +
  scale_colour_manual(values = method_col, name = NULL) +
  labs(x = "Mean SE / RMSE", y = NULL, tag = "b") +
  theme_paper() + theme(strip.placement = "outside", strip.text = element_blank(), axis.text.y = element_blank(),
                        panel.grid.major.y = element_blank()) +
  legend_inside(0.98, 0.02, c(1, 0))

## 2. Panel c: coverage by mixture design ---------------------------------------------------
# A second simulation whose mixtures copy the separations of the empirical fits; each
# design has several cells of sample size, allele frequency and effect.
lowsep <- as.data.table(read_result("08_simulations/summary", "lowsep_sim_summary"))
geom_lab <- c(testosterone_unshrunk = "Testosterone,\nunshrunk", testosterone_shrink_0.2 = "Shrink 0.2", testosterone_shrink_0.5 = "Shrink 0.5",
              testosterone_shrink_0.7 = "Shrink 0.7", mcv_vp1 = "MCV VP1", jcv_vp1 = "JCV VP1")
lowsep <- lowsep[geometry %in% names(geom_lab)]
cov_cells <- rbind(
  lowsep[, .(geometry, separation, parameter, n, maf, effect, cov = "G-MoQ, pooled", coverage = pooled_coverage95, mcse = pooled_coverage95_mcse)],
  lowsep[, .(geometry, separation, parameter, n, maf, effect, cov = "G-MoQ, jackknife", coverage = jackknife_coverage95, mcse = jackknife_coverage95_mcse)])
cov_sum <- cov_cells[, .(coverage = mean(coverage), lo = min(coverage), hi = max(coverage), cells = .N), by = .(geometry, parameter, cov)]
sep <- unique(lowsep[, .(geometry, separation)])
cov_sum <- merge(cov_sum, sep, by = "geometry")
geom_row <- function(g, s) sprintf("%s (S %.1f)", sub("\n", " ", geom_lab[g]), s)
cov_sum[, g_f := factor(geom_row(geometry, separation), levels = sep[order(-separation), geom_row(geometry, separation)])]
cov_sum[, par_f := factor(c(beta_1 = "delta[1]", beta_2 = "delta[2]", gamma = "gamma")[parameter], levels = unname(par_lab))]
cov_sum[, cov := factor(cov, levels = method_lab[1:2])]
median_mcse <- median(cov_cells$mcse)
pC <- ggplot(cov_sum, aes(coverage, g_f, colour = cov)) +
  annotate("rect", xmin = 0.95 - 1.96 * median_mcse, xmax = 0.95 + 1.96 * median_mcse, ymin = -Inf, ymax = Inf, fill = "grey90") +
  geom_vline(xintercept = 0.95, linewidth = 0.3, linetype = 2, colour = "grey45") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, linewidth = 0.5, position = position_dodge(width = 0.5), orientation = "y") +
  geom_point(size = 1.6, position = position_dodge(width = 0.5)) +
  facet_wrap(~par_f, nrow = 1, labeller = label_parsed) +
  scale_colour_manual(values = method_col, guide = "none") +
  scale_x_continuous(breaks = c(0.92, 0.95, 0.98)) +
  labs(x = "95% CI coverage", y = NULL, tag = "c") +
  theme_paper() + theme(panel.grid.major.y = element_blank())

## 3. Panel d: time per 100-variant panel ---------------------------------------------------
timings <- as.data.table(read_result("08_simulations/summary", "mixture_timings"))
timings[, method_f := factor(method_lab[method], levels = method_lab)]
pD <- ggplot(timings, aes(n, median_fit_seconds, colour = method_f)) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.5) +
  scale_x_log10(breaks = c(1e4, 5e4, 4e5), labels = c("10k", "50k", "400k")) +
  scale_y_log10(breaks = c(3, 10, 30, 100, 300, 1000), labels = c("3", "10", "30", "100", "300", "1,000")) +
  scale_colour_manual(values = method_col, guide = "none") +
  labs(x = "Sample size", y = "Fit + SE time per panel (s)", tag = "d") +
  theme_paper()

## 4. Assemble and save ---------------------------------------------------------------------
fig <- ((pA | pB) + plot_layout(widths = c(1.3, 1))) / free((pC | pD) + plot_layout(widths = c(2.4, 1))) + plot_layout(heights = c(1.5, 1))
src <- rbind(sim[, .(panel = "a/b", design = "S4", parameter, n, method, value1 = bias, value2 = mean_se / rmse, value3 = NA_real_)],
             cov_sum[, .(panel = "c", design = geometry, parameter, n = NA_real_, method = as.character(cov), value1 = coverage, value2 = lo, value3 = hi)],
             timings[, .(panel = "d", design = "M100", parameter = NA_character_, n, method, value1 = median_fit_seconds, value2 = NA_real_, value3 = NA_real_)])
t_at <- function(nn, m) timings[n == nn & method == m, median_fit_seconds]
n_pan <- unique(sim$success)
cells_per_design <- lowsep[, .(k = uniqueN(cell)), by = geometry]
cap <- paste(
  "Figure S4. G-MoQ Stage 2 and direct mixture maximum likelihood on simulated two-component mixtures.",
  "Parameters are the lower and upper component locations (δ_1, δ_2) and membership (γ), each against the",
  "method's own known truth; G-MoQ with pooled or 1,000-block jackknife covariance shares point estimates (Methods).",
  sprintf(paste("Panels a and b use the same simulated datasets at three sample sizes: %s replicate panels of 100",
                "variants per sample size. **a**, Bias; point, mean over the panels; error bars, 95%% Monte Carlo",
                "intervals (±1.96 RMSE/√%s). **b**, Mean SE over the RMSE; the dashed line marks 1."), n_pan, n_pan),
  sprintf(paste("**c**, Coverage of the nominal 95%% confidence interval in a second simulation, by mixture design",
                "(weights fitted to the testosterone phenotype, unshrunk and shrunk, to JCV VP1 and to MCV VP1;",
                "S, component separation). Each design has %s cells of sample size, minor allele frequency",
                "and effect, each of %s replicates; point, mean over cells; bar, range over cells; grey band,",
                "0.95 ± 1.96 times the median Monte Carlo standard error; dashed line, 0.95."),
          paste(sort(unique(cells_per_design$k)), collapse = " or "), fmt(unique(lowsep$replicates_ok))),
  sprintf(paste("**d**, Median time to fit one 100-variant panel with standard errors, over %d panels per sample size:",
                "G-MoQ pooled %.1f s against mixture maximum likelihood %.0f s at 10,000 participants, and",
                "%.0f s against %s s at 400,000."), unique(timings$n_panels), t_at(1e4, "gmoq_pooled"), t_at(1e4, "mixture_ml"),
          t_at(4e5, "gmoq_pooled"), fmt(round(t_at(4e5, "mixture_ml")))),
  "No statistical tests are shown.",
  "ML, maximum likelihood; RMSE, root mean square error; SE, standard error; CI, confidence interval;",
  "JCV, JC polyomavirus; MCV, Merkel cell polyomavirus; VP1, major capsid protein VP1.")
save_figure(fig, "figS04_mixture_recovery", src, cap, width = 183, height = 150,
            results_used = "08_simulations/summary")
