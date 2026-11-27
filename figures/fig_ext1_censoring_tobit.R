#!/usr/bin/env Rscript
# Extended Data / Supplementary Figure: Limit-of-detection censoring via G-MoQ vs Tobit MLE.
#   a  Simulation benchmark (200 replicates, n = 20,000) across three censoring regimes:
#      Gaussian double-censored (11% L, 7% R), Gaussian heavy left-censored (70% L),
#      and bimodal double-censored (11% L, 7% R). Relative bias (left) and mean |z| (right).
#   b  Real UK Biobank female testosterone (97 female lead variants; 15.7% left-censored
#      below 0.35 nmol/L): estimated log-testosterone effect vs returned fitted Tobit estimates.
suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(data.table)})
.repo <- system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE)
source(file.path(.repo, "figures", "theme.R"))

# ---- Panel a: Simulation benchmark ----------------------------------------------------------
sim <- fread(file.path(.repo, "results", "03_testosterone", "censoring_tobit", "lod_tobit_sim.csv"))
scen_map <- c(
  S1_gauss_lpa_like_11L_7R = "Gaussian\n(11% L, 7% R)",
  S2_gauss_uacr_like_70L   = "Gaussian\n(70% L)",
  S3_bimodal_lpa_like_11L_7R = "Bimodal\n(11% L, 7% R)"
)
meth_map <- c(
  drop_ols  = "Drop OLS",
  clamp_ols = "Clamp OLS",
  tobit     = "Tobit MLE",
  gmoq_ols  = "Censored G-MoQ (OLS)"
)
# OLS with pooled covariance is the primary estimator.
mcol_a <- c(
  "Drop OLS"             = "#9aa3b2",
  "Clamp OLS"            = "#c8a26a",
  "Tobit MLE"      = "#3d4c63",
  "Censored G-MoQ (OLS)" = paper_palette[["shift"]]
)

sim_p <- sim[method %in% names(meth_map)]
sim_p[, scen_f := factor(scen_map[scenario], levels = unname(scen_map))]
sim_p[, meth_f := factor(meth_map[method], levels = unname(meth_map))]
orc <- sim[method == "oracle"]
orc[, scen_f := factor(scen_map[scenario], levels = unname(scen_map))]

pA1 <- ggplot(sim_p, aes(scen_f, rel_bias_pct, fill = meth_f)) +
  geom_hline(yintercept = 0, linewidth = 0.35, colour = "grey35") +
  geom_col(position = position_dodge(width = 0.78), width = 0.70) +
  scale_fill_manual(values = mcol_a, name = NULL) +
  scale_y_continuous(limits = c(-85, 24), breaks = c(-80, -60, -40, -20, 0, 20)) +
  labs(x = NULL, y = "Relative bias (%)", tag = "a") +
  theme_paper() +
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.key.size = unit(0.55, "lines"), legend.text = element_text(size = 5.5),
        legend.margin = margin(0, 0, 2, 0)) +
  guides(fill = guide_legend(nrow = 2))

pA2 <- ggplot(sim_p, aes(scen_f, mean_abs_z, fill = meth_f)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.70) +
  geom_errorbar(data = orc, aes(x = scen_f, ymin = mean_abs_z, ymax = mean_abs_z),
                inherit.aes = FALSE, width = 0.80, linetype = "dashed", linewidth = 0.55, colour = "grey15") +
  scale_fill_manual(values = mcol_a, name = NULL, guide = "none") +
  scale_y_continuous(limits = c(0, 4.7), breaks = 0:4, expand = expansion(mult = c(0, 0.03))) +
  labs(x = NULL, y = "Mean |z|") +
  theme_paper()

# ---- Panel b: UK Biobank female testosterone (97 leads) ------------------------------------
tes <- fread(file.path(.repo, "results", "03_testosterone", "censoring_tobit", "testo_tobit_gmoq_results.csv"))[is_fem == TRUE]
# Orient each SNP so its fitted Tobit effect is positive for display.
# Agreement statistics in the caption use the original signed estimates.
sgn <- ifelse(tes$b_tobit >= 0, 1, -1)
tb <- rbind(
  data.table(snp_id = tes$snp_id, display_sign = sgn, signed_x = tes$b_tobit,
             x = tes$b_tobit * sgn, y = tes$b_f_cens_ols * sgn, signed_y = tes$b_f_cens_ols,
             method = "Censored G-MoQ OLS (female)"),
  data.table(snp_id = tes$snp_id, display_sign = sgn, signed_x = tes$b_tobit,
             x = tes$b_tobit * sgn, y = tes$b_mix_cens_f * sgn, signed_y = tes$b_mix_cens_f,
             method = "Sex-blind censored mixture OLS"),
  data.table(snp_id = tes$snp_id, display_sign = sgn, signed_x = tes$b_tobit,
             x = tes$b_tobit * sgn, y = tes$b_f_drop * sgn, signed_y = tes$b_f_drop,
             method = "Dropped OLS (female)"),
  data.table(snp_id = tes$snp_id, display_sign = sgn, signed_x = tes$b_tobit,
             x = tes$b_tobit * sgn, y = tes$b_mix_drop_f * sgn, signed_y = tes$b_mix_drop_f,
             method = "Sex-blind dropped mixture OLS")
)
blevels <- c("Censored G-MoQ OLS (female)", "Sex-blind censored mixture OLS",
             "Dropped OLS (female)", "Sex-blind dropped mixture OLS")
tb[, meth_f := factor(method, levels = blevels)]
bcol <- c(
  "Censored G-MoQ OLS (female)"    = paper_palette[["shift"]],
  "Sex-blind censored mixture OLS" = paper_palette[["spread"]],
  "Dropped OLS (female)"           = "#8d96a5",
  "Sex-blind dropped mixture OLS"  = "#c8a26a"
)
bshp <- c(
  "Censored G-MoQ OLS (female)"    = 16,
  "Sex-blind censored mixture OLS" = 17,
  "Dropped OLS (female)"           = 1,
  "Sex-blind dropped mixture OLS"  = 2
)

pB <- ggplot(tb, aes(x, y, colour = meth_f, shape = meth_f)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.45, colour = "grey30") +
  geom_point(size = 1.25, alpha = 0.80, stroke = 0.45) +
  scale_colour_manual(values = bcol, name = NULL) +
  scale_shape_manual(values = bshp, name = NULL) +
  coord_cartesian(xlim = c(0, 0.135), ylim = c(0, 0.135)) +
  scale_x_continuous(breaks = c(0, 0.04, 0.08, 0.12)) +
  scale_y_continuous(breaks = c(0, 0.04, 0.08, 0.12)) +
  labs(x = "Fitted female Tobit effect (log nmol/L)",
       y = "Estimated female effect (log nmol/L)", tag = "b") +
  theme_paper() + legend_inside(0.02, 0.98, c(0, 1)) +
  theme(legend.key.size = unit(0.65, "lines"))

# ---- Combine panels a (sim bias + mean |z|) and b (UKB female testosterone) -----------------
p <- (pA1 | pA2 | pB) + plot_layout(widths = c(1.05, 0.95, 1.15))

src <- rbindlist(list(
  sim_p[, .(panel = "a", series = paste(scenario, method), x = rel_bias_pct, y = mean_abs_z)],
  tb[, .(panel = "b", series = method, snp_id, x, y, display_sign, signed_x, signed_y)]
), fill = TRUE)
cap <- paste(
  "Extended Data Figure 1. Limit-of-detection censoring in G-MoQ and fitted Tobit benchmarks.",
  "a, Simulation benchmark (200 replicates per scenario, n = 20,000, true beta = 0.04) comparing relative bias (left)",
  "and mean absolute z statistic (right) across double-censored Gaussian (11% left, 7% right), heavily left-censored",
  "Gaussian (70% left), and double-censored bimodal (11% left, 7% right) distributions; dashed horizontal bars mark",
  "uncensored oracle ordinary least squares. b, Real UK Biobank female testosterone lead variants (n = 97; 199,546",
  "observed and 37,037 left-censored women below 0.35 nmol/L, 15.7% censored), comparing sex-stratified and sex-blind",
  "G-MoQ estimators against the returned female Tobit estimates (dashed identity line; censored G-MoQ slope = 0.992,",
  "Pearson r = 0.9996; dropped OLS slope = 0.578; sex-blind censored mixture slope = 1.002, r = 0.995). Agreement statistics use the original signed estimates; each SNP is oriented to a positive Tobit effect for display. A successful optimiser exit was recorded at seven of the 97 leads.")

save_figure(p, "fig_ext1_censoring_tobit", src, cap, width = 183, height = 78,
            results_used = c("03_testosterone/censoring_tobit"))
