#!/usr/bin/env Rscript
# Fig S4: G-MoQ Stage 2 (pooled, jackknife) vs full mixture maximum likelihood on the same simulated
# datasets. Replaces the beta / SE scatter with bias, SE / error SD and coverage by parameter and
# design (author review, 28 Sep 2026). Reads results/P17_simulations only:
#   S4 rows of table_s2_numerical_summary (bias, rmse, mean_se; 20 replicate panels per N; no coverage
#   recorded), lowsep_sim_summary (A2 mixture simulation: coverage by parameter and design) and
#   mixture_timings.
suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(data.table)})
.repo <- system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE)
source(file.path(.repo, "figures", "theme.R"))
x <- as.data.table(read_result("08_simulations/summary", "table_s2_numerical_summary"))[study == "S4"]
mlab <- c(gmoq_pooled = "G-MoQ, pooled", gmoq_jackknife = "G-MoQ, jackknife", mixture_ml = "Mixture ML")
mcol <- setNames(c(paper_palette[["gmoq"]], "#7fb3d5", paper_palette[["mixture"]]), mlab)
x[, method_f := factor(mlab[method], levels = mlab)]
x[, n_f := factor(sprintf("N = %s", format(n, big.mark = ",", scientific = FALSE, trim = TRUE)),
                  levels = sprintf("N = %s", format(sort(unique(n)), big.mark = ",", scientific = FALSE, trim = TRUE)))]
# Greek as plotmath (the Mac PDF fallback device cannot draw Unicode Greek)
plab <- c(delta1 = "delta[1]", delta2 = "delta[2]", gamma = "gamma")
x[, par_f := factor(plab[parameter], levels = plab)]
x[, mcse := rmse / sqrt(success)]
pd <- position_dodge(width = 0.6)
pA <- ggplot(x, aes(bias, n_f, colour = method_f)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_errorbar(aes(xmin = bias - 1.96 * mcse, xmax = bias + 1.96 * mcse), width = 0, linewidth = 0.6, position = pd,
                orientation = "y") +
  geom_point(size = 1.7, position = pd) +
  facet_wrap(~par_f, ncol = 1, strip.position = "left", labeller = label_parsed) +
  scale_colour_manual(values = mcol, name = NULL) +
  labs(x = "Bias", y = NULL, tag = "a") +
  theme_paper() + theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1, size = rel(1.15)),
                        panel.grid.major.y = element_blank(), legend.position = "none")
pB <- ggplot(x, aes(mean_se / rmse, n_f, colour = method_f)) +
  geom_vline(xintercept = 1, linewidth = 0.3, linetype = 2, colour = "grey45") +
  geom_point(size = 1.7, position = pd) +
  facet_wrap(~par_f, ncol = 1, strip.position = "left", labeller = label_parsed) +
  scale_colour_manual(values = mcol, name = NULL) +
  labs(x = "Mean SE / RMSE", y = NULL, tag = "b") +
  theme_paper() + theme(strip.placement = "outside", strip.text = element_blank(), axis.text.y = element_blank(),
                        panel.grid.major.y = element_blank()) +
  legend_inside(0.98, 0.02, c(1, 0))
# coverage by parameter and mixture design (A2 low-separation simulation)
ls <- as.data.table(read_result("08_simulations/summary", "lowsep_sim_summary"))
glab <- c(testosterone_unshrunk = "Testosterone,\nunshrunk", testosterone_shrink_0.2 = "Shrink 0.2", testosterone_shrink_0.5 = "Shrink 0.5",
          testosterone_shrink_0.7 = "Shrink 0.7", mcv_vp1 = "MCV VP1", jcv_vp1 = "JCV VP1")
ls <- ls[geometry %in% names(glab)]
cvl <- rbind(ls[, .(geometry, separation, parameter, n, maf, effect, cov = "G-MoQ, pooled", coverage = pooled_coverage95, mcse = pooled_coverage95_mcse)],
             ls[, .(geometry, separation, parameter, n, maf, effect, cov = "G-MoQ, jackknife", coverage = jackknife_coverage95, mcse = jackknife_coverage95_mcse)])
cvs <- cvl[, .(coverage = mean(coverage), lo = min(coverage), hi = max(coverage), cells = .N), by = .(geometry, parameter, cov)]
sep <- unique(ls[, .(geometry, separation)])
cvs <- merge(cvs, sep, by = "geometry")
cvs[, g_f := factor(sprintf("%s (S %.1f)", sub("\n", " ", glab[geometry]), separation),
                    levels = sep[order(-separation), sprintf("%s (S %.1f)", sub("\n", " ", glab[geometry]), separation)])]
cvs[, par_f := factor(c(beta_1 = "delta[1]", beta_2 = "delta[2]", gamma = "gamma")[parameter], levels = unname(plab))]
cvs[, cov := factor(cov, levels = mlab[1:2])]
mc0 <- median(cvl$mcse)
pC <- ggplot(cvs, aes(coverage, g_f, colour = cov)) +
  annotate("rect", xmin = 0.95 - 1.96 * mc0, xmax = 0.95 + 1.96 * mc0, ymin = -Inf, ymax = Inf, fill = "grey90") +
  geom_vline(xintercept = 0.95, linewidth = 0.3, linetype = 2, colour = "grey45") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, linewidth = 0.5, position = position_dodge(width = 0.5), orientation = "y") +
  geom_point(size = 1.6, position = position_dodge(width = 0.5)) +
  facet_wrap(~par_f, nrow = 1, labeller = label_parsed) +
  scale_colour_manual(values = mcol, guide = "none") +
  scale_x_continuous(breaks = c(0.92, 0.95, 0.98)) +
  labs(x = "95% CI coverage", y = NULL, tag = "c") +
  theme_paper() + theme(panel.grid.major.y = element_blank())
tm_raw <- as.data.table(read_result("08_simulations/speed_benchmark", "mixture_panel_repeats"))[method %in% names(mlab)]
tm <- tm_raw[, .(median_fit_seconds = median(fit_seconds), n_panels = .N), by = .(n, method)]
tm[, method_f := factor(mlab[method], levels = mlab)]
pD <- ggplot(tm, aes(n, median_fit_seconds, colour = method_f)) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.5) +
  scale_x_log10(breaks = c(1e4, 5e4, 4e5), labels = c("10k", "50k", "400k")) +
  scale_y_log10(breaks = c(1, 3, 10, 30, 100, 300, 1000), labels = c("1", "3", "10", "30", "100", "300", "1,000")) +
  scale_colour_manual(values = mcol, guide = "none") +
  labs(x = "Sample size", y = "Fit + SE time per panel (s)", tag = "d") +
  theme_paper()
p <- ((pA | pB) + plot_layout(widths = c(1.3, 1))) / free((pC | pD) + plot_layout(widths = c(2.4, 1))) + plot_layout(heights = c(1.5, 1))
src <- rbind(x[, .(panel = "a/b", design = "S4", parameter, n, method, value1 = bias, value2 = mean_se / rmse, value3 = NA_real_)],
             cvs[, .(panel = "c", design = geometry, parameter, n = NA_real_, method = as.character(cov), value1 = coverage, value2 = lo, value3 = hi)],
             tm[, .(panel = "d", design = "M100", parameter = NA_character_, n, method, value1 = median_fit_seconds, value2 = NA_real_, value3 = NA_real_)])
r10 <- tm[n == 1e4 & method == "mixture_ml", median_fit_seconds]; g10 <- tm[n == 1e4 & method == "gmoq_pooled", median_fit_seconds]
r400 <- tm[n == 4e5 & method == "mixture_ml", median_fit_seconds]; g400 <- tm[n == 4e5 & method == "gmoq_pooled", median_fit_seconds]
fmt <- function(v) format(v, big.mark = ",", scientific = FALSE, trim = TRUE)
n_pan <- unique(x$success); n_ls <- unique(ls$replicates_ok)
ncell <- ls[, .(k = uniqueN(cell)), by = geometry]
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
          paste(sort(unique(ncell$k)), collapse = " or "), fmt(n_ls)),
  sprintf(paste("**d**, Median time to fit one 100-variant panel with standard errors, over %d panels per sample size:",
                "G-MoQ pooled %.1f s against mixture maximum likelihood %.0f s at 10,000 participants, and",
                "%.0f s against %s s at 400,000."), unique(tm$n_panels), g10, r10, g400, fmt(round(r400))),
  "No statistical tests are shown.",
  "ML, maximum likelihood; RMSE, root mean square error; SE, standard error; CI, confidence interval;",
  "JCV, JC polyomavirus; MCV, Merkel cell polyomavirus; VP1, major capsid protein VP1.")
save_figure(p, "figS04_mixture_recovery", src, cap, width = 183, height = 150,
            results_used = c("08_simulations/summary/table_s2_numerical_summary", "08_simulations/summary/lowsep_sim_summary",
                             "08_simulations/summary/mixture_timings"))
