#!/usr/bin/env Rscript
# Figure 1: reading a sex-blind quantile GWAS of testosterone through a two-component mixture
#
# What: testosterone is a mixture of a female and a male distribution. Stage 1 sees only the
#   whole-cohort distribution; Stage 2 maps its quantile effects onto the parameters of a
#   fitted two-component skew-normal mixture through the weights W[tau, p]. The figure shows
#   the distributions being decomposed and the weights that do it, for the unshrunk trait and
#   three versions in which the male values are shrunk towards the female mean.
#   a: distribution of the analysed phenotype, whole cohort and its female and male parts.
#   b: Stage 2 weights for female location, male location and membership (log odds of the upper component).
# Paper: Fig. 1
# Input: results/03_testosterone/recovery/phenotype_histograms.csv,
#   results/02_mixtures/stage2_mixtures/weight_matrices.csv
# Output: paper_assets/fig1_testosterone_mixture.{png,pdf}, _source.csv, _caption.md
# Runs: locally from results/ in about 2 seconds

suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(data.table)})
source(file.path(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE), "figures", "theme.R"))

cfg <- read_config()
arm_lab <- c(testosterone = "No shrinkage", testosterone_shrink_0_2 = "Shrink 0.2",
             testosterone_shrink_0_5 = "Shrink 0.5", testosterone_shrink_0_7 = "Shrink 0.7")
tau_axis <- scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", ".25", ".5", ".75", "1"))

## 1. Load ----------------------------------------------------------------------------------
hist <- as.data.table(read_result("03_testosterone/recovery", "phenotype_histograms"))
w_all <- as.data.table(read_result("02_mixtures/stage2_mixtures", "weight_matrices"))[fit %in% names(arm_lab)]

## 2. Panel a: phenotype distributions -------------------------------------------------------
hist[, arm := gsub(".", "_", phenotype, fixed = TRUE)]
hist <- hist[is.finite(bin_lower) & is.finite(bin_upper) & !is.na(density)]
n_all <- hist[group == "all", .(n_group = n_group[1]), by = arm][, setNames(n_group, arm)]
# Each sex's histogram is scaled by that sex's share of the cohort, so female + male stack
# to the whole-cohort distribution.
hist[, dens_plot := density * ifelse(group == "all", 1, n_group / n_all[arm])]
hist[, arm_f := factor(arm_lab[arm], levels = arm_lab)]
# Display densities: each histogram convolved with a Gaussian kernel of one median bin width.
grid <- seq(-1.1, 3.4, length.out = 500)
dens <- hist[, { bw <- 1.0 * median(bin_upper - bin_lower)
                 .(x = grid, y = smooth_hist_density(bin_lower, bin_upper, dens_plot, grid, bw)) },
             by = .(arm, arm_f, group)]
dens <- dens[y > 2e-3]
pA <- ggplot() +
  geom_area(data = dens[group != "all"], aes(x, y, fill = group), alpha = 0.5, position = "identity",
            colour = NA, outline.type = "full") +
  geom_line(data = dens[group != "all"], aes(x, y, colour = group), linewidth = 0.35) +
  geom_line(data = dens[group == "all"], aes(x, y, colour = group), linewidth = 0.6, lineend = "round") +
  facet_wrap(~arm_f, nrow = 1) +
  scale_fill_manual(values = c(female = paper_palette[["female"]], male = paper_palette[["male"]]),
                    labels = c(female = "Female", male = "Male"), name = NULL) +
  scale_colour_manual(values = c(all = "grey10", female = paper_palette[["female"]], male = paper_palette[["male"]]),
                      breaks = "all", labels = c(all = "Whole cohort"), name = NULL) +
  scale_x_continuous(breaks = -1:3) +
  guides(fill = guide_legend(order = 1), colour = guide_legend(order = 2)) +
  labs(x = "Testosterone (analysed scale)", y = "Density", tag = "a") +
  theme_paper() + legend_inside(0.02, 0.98, c(0, 1)) + theme(legend.key.width = unit(0.9, "lines"))

## 3. Panel b: Stage 2 weights ---------------------------------------------------------------
# Each arm is drawn with its reported weight design (config.yaml): 'hybrid' for the unshrunk arm
# and F = 0.2, 'model' otherwise. Where the hybrid design is used, the model-design membership
# weight is drawn as well (solid, like the other arms) and the hybrid weight dashed.
# Membership is shown as the log odds of the UPPER (male) component, i.e. -gamma: the code's
# gamma is the log odds of the lower component.
w <- pick_design(w_all, "fit", cfg = cfg)
w[, arm_f := factor(arm_lab[fit], levels = arm_lab)]
hybrid_fits <- unique(w[weight_design == "hybrid"]$fit)
w_model_gamma <- w_all[fit %in% hybrid_fits & weight_design == "model" & param == "gamma"]
w_model_gamma[, arm_f := factor(arm_lab[fit], levels = arm_lab)]
# The weights depend only on the baseline mixture, so they are the same on every chromosome;
# chromosome 1 is drawn.
stopifnot(max(abs(w$W_max_chr - w$W_min_chr)) < 1e-6)

w_lines <- w[!(fit %in% hybrid_fits & param == "gamma")]
w_gamma_both <- rbind(w_model_gamma[, .(arm_f, param, tau, W_chr01, line = "model")],
                      w[fit %in% hybrid_fits & param == "gamma", .(arm_f, param, tau, W_chr01, line = "hybrid")])
curve_df <- function(d, extra = character(0)) smooth_curve(d, "tau", "W_chr01", by = c("arm_f", "param", extra), n = 600)
loc_curves <- curve_df(w_lines[param != "gamma"])
gamma_curves <- curve_df(w_lines[param == "gamma"]); gamma_curves[, line := "model"]
gamma_curves <- rbind(gamma_curves, curve_df(w_gamma_both, "line"))
gamma_curves[, y := -y]                        # membership of the upper component
row_strips <- theme(strip.text.y.right = element_text(angle = 0, hjust = 0, face = "bold"))

# Top rows: female and male location share one 0-1 scale.
loc_curves[, param_f := factor(param, levels = c("beta_left", "beta_right"), labels = c("Female\nlocation", "Male\nlocation"))]
pB_loc <- ggplot(loc_curves, aes(x, y, colour = param)) +
  geom_hline(yintercept = c(0, 1), linewidth = 0.25, linetype = 2, colour = "grey60") +
  geom_line(linewidth = 0.7, lineend = "round") +
  facet_grid(param_f ~ arm_f) + tau_axis +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  scale_colour_manual(values = c(beta_left = paper_palette[["female"]], beta_right = paper_palette[["male"]]), guide = "none") +
  labs(x = NULL, y = "Stage 2 weight", tag = "b") +
  theme_paper() + row_strips + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
# Bottom row: membership on its own scale.
gamma_curves[, param_f := "Membership\n(gamma)"]
pB_gamma <- ggplot(gamma_curves, aes(x, y)) +
  geom_hline(yintercept = 0, linewidth = 0.25, linetype = 2, colour = "grey60") +
  geom_line(data = gamma_curves[line == "model"], colour = paper_palette[["mixture"]], linewidth = 0.7, lineend = "round") +
  geom_line(data = gamma_curves[line == "hybrid"], colour = "grey10", linewidth = 0.6, linetype = "22") +
  facet_grid(param_f ~ arm_f) + tau_axis +
  labs(x = "Quantile level tau", y = NULL) +
  theme_paper() + row_strips + theme(strip.text.x = element_blank())
pB <- pB_loc / pB_gamma + plot_layout(heights = c(2, 1.15))

## 4. Assemble and save ---------------------------------------------------------------------
fig <- (pA / pB) + plot_layout(heights = c(1, 1.9))
src <- rbind(
  hist[group != "all", .(panel = "A", arm, series = group, x_lower = bin_lower, x_upper = bin_upper, value = dens_plot, n = n_group)],
  hist[group == "all", .(panel = "A", arm, series = "all", x_lower = bin_lower, x_upper = bin_upper, value = dens_plot, n = n_group)],
  w[, .(panel = "B", arm = fit, series = paste(param, weight_design), x_lower = tau, x_upper = tau,
        value = ifelse(param == "gamma", -W_chr01, W_chr01), n = NA_real_)],
  w_model_gamma[, .(panel = "B", arm = fit, series = paste(param, weight_design), x_lower = tau, x_upper = tau, value = -W_chr01, n = NA_real_)])
cap <- paste(
  "Figure 1. Reading a sex-blind quantile GWAS through a two-component mixture model.",
  "Testosterone, UK Biobank; n = 399,710 (199,572 female, 200,138 male).",
  "**a**, Distribution of the analysed phenotype in the unshrunk trait and in three arms in which the male-female",
  "separation is shrunk by 0.2, 0.5 and 0.7, increasing the overlap. The black line is the whole-cohort distribution",
  "seen by the sex-blind Stage 1 analysis; shaded areas are its female and male parts, each scaled by its share.",
  "**b**, Stage 2 weights, the derivative of each quantile with respect to each parameter, for the fitted",
  "two-component skew-normal mixture: female location, male location and membership gamma, the log odds of",
  "belonging to the upper (male) component. Solid lines use the model-density weights; dashed lines, in the unshrunk",
  "arm and at 0.2, are the hybrid weights used there, whose membership column is divided by the Stage 1 kernel density.",
  "Weights are shown for chromosome 1 and are identical on every chromosome. Stage 2 is ordinary least squares with",
  "the pooled Stage 1 covariance (Methods).")
save_figure(fig, "fig1_testosterone_mixture", src, cap, width = 200, height = 150,
            results_used = c("03_testosterone/recovery", "02_mixtures/stage2_mixtures"))
