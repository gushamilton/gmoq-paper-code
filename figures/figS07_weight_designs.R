#!/usr/bin/env Rscript
# Supplementary Figure S7: the model and hybrid Stage 2 weight designs for two-component mixtures
#
# What: Stage 2 mixture weights divide by a density. The model design uses the fitted mixture
#   density for every column; the hybrid design uses the Stage 1 kernel density estimate (KDE)
#   for the membership (gamma) column. Where the fitted mixture misses the density in the gap
#   between components (unshrunk testosterone), the two differ. Panel a shows where; panel b
#   shows, by simulation, how each design recovers a membership-only effect.
#   a: fitted mixture density over the Stage 1 KDE across quantile levels, testosterone arms.
#   b: mean gamma estimate over its true value in nine scenarios, with coverage per design.
# Paper: Supplementary Fig. S7
# Input: results/02_mixtures/stage2_mixtures/weight_density.csv,
#   results/02_mixtures/weight_design/sim/{operating_characteristics,fit_summary}.csv
# Output: paper_assets/figS07_weight_designs.{png,pdf}, _source.csv, _caption.md
# Runs: locally from results/ in about 2 seconds

suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(data.table)})
source(file.path(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE), "figures", "theme.R"))

arm_shrink <- c(testosterone = 0, testosterone_shrink_0_2 = 0.2, testosterone_shrink_0_5 = 0.5,
                testosterone_shrink_0_7 = 0.7)
arm_lab <- c("0", "0.2", "0.5", "0.7")
arm_col <- c("#3b1f5c", "#3f6ea8", "#2a9d8f", "#9bbf3a")
design_col <- c(hybrid = paper_palette[["gmoq"]], model = paper_palette[["mixture"]])
design_lab <- c(hybrid = "Hybrid", model = "Model")

## 1. Panel a: mixture density over the Stage 1 KDE ----------------------------------------
wd <- as.data.table(read_result("02_mixtures/stage2_mixtures", "weight_density"))
wd <- wd[fit %in% names(arm_shrink)]
wd[, arm_f := factor(arm_lab[match(fit, names(arm_shrink))], levels = arm_lab)]
gap <- wd[abs(tau - 0.5) < 1e-6]                 # tau = 0.5 lies between the two modes
stopifnot(nrow(gap) == 4)
n_tau <- wd[fit == "testosterone", .N]
pA <- ggplot(wd, aes(tau, ratio_model_over_kde, colour = arm_f)) +
  geom_hline(yintercept = 1, linetype = 2, linewidth = 0.25, colour = "grey50") +
  geom_line(data = smooth_curve(wd, "tau", "ratio_model_over_kde", by = "arm_f", n = 600), aes(x, y), linewidth = 0.65) +
  geom_point(data = gap, size = 1.5) +
  scale_colour_manual(values = setNames(arm_col, arm_lab), name = "Shrinkage") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", ".25", ".5", ".75", "1")) +
  labs(x = expression("Quantile level " * tau), y = "Mixture density / Stage 1 KDE", tag = "a") +
  theme_paper() + legend_inside(0.02, 0.03, c(0, 0)) +
  theme(legend.title = element_text(size = 7))

## 2. Panel b: simulated recovery of a membership-only effect -----------------------------
# Design codes in the simulation output: 'par' = model design, 'hyb' = hybrid design.
oc <- as.data.table(read_result("02_mixtures/weight_design/sim", "operating_characteristics"))
fs <- as.data.table(read_result("02_mixtures/weight_design/sim", "fit_summary"))
gam <- oc[role == "gamma" & param == "gamma" & design %in% c("par", "hyb")]
gam[, `:=`(ratio = mean_est / truth, ratio_se = mcse_bias / truth,
           weight_design = ifelse(design == "hyb", "hybrid", "model"))]
scen_lab <- c(ctrl_sep = "separated, correct mixture", ctrl_sep_n400k = "separated, correct, n = 400k",
              heavy_c06 = "testosterone-like gap misfit", heavy_c06_n400k = "testosterone-like, n = 400k",
              heavy_c03 = "heavier gap misfit", heavy_c06_p03 = "testosterone-like, lower share 0.3",
              bimodal_sep = "bimodal separated", heavy_c06_shrink05 = "overlap (shrink 0.5)",
              heavy_c06_shrink07 = "heavy overlap (shrink 0.7)")
gam <- gam[scenario %in% names(scen_lab)]
stopifnot(nrow(gam) == 2 * length(scen_lab))
gam[, sc_f := factor(scen_lab[scenario], levels = rev(scen_lab))]
coverage <- dcast(gam, sc_f ~ weight_design, value.var = "coverage")
pB <- ggplot(gam, aes(ratio, sc_f, colour = weight_design)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.25, colour = "grey50") +
  geom_errorbar(aes(xmin = ratio - 1.96 * ratio_se, xmax = ratio + 1.96 * ratio_se), width = 0.3, orientation = "y") +
  geom_point(size = 1.4) +
  # Coverage of nominal 95% intervals, printed as a labelled column right of the panel.
  geom_text(data = coverage, aes(x = 1.07, y = sc_f, label = sprintf("%.2f / %.2f", model, hybrid)), inherit.aes = FALSE,
            hjust = 0, size = 6 / .pt, colour = "grey20") +
  annotate("text", x = 1.07, y = length(scen_lab) + 0.75, label = "Coverage\nmodel / hybrid", hjust = 0,
           size = 6 / .pt, lineheight = 0.9, colour = "grey20") +
  scale_colour_manual(values = design_col, labels = design_lab, name = NULL) +
  scale_x_continuous(limits = c(0.3, 1.32), breaks = seq(0.4, 1, 0.2)) +
  coord_cartesian(clip = "off") +
  labs(x = expression(gamma * " estimate / true " * gamma), y = NULL, tag = "b") +
  theme_paper() + legend_inside(0.02, 0.02, c(0, 0)) +
  theme(plot.margin = margin(14, 6, 4, 4))

## 3. Numbers for the legend ----------------------------------------------------------------
# Location bias on the narrower (upper) component for a variant acting on it alone, in the
# seven separated scenarios, by design.
fs <- fs[scenario %in% names(scen_lab)]
stopifnot(all(fs$omega2 < fs$omega1))            # the upper component is the narrower one
separated <- setdiff(names(scen_lab), c("heavy_c06_shrink05", "heavy_c06_shrink07"))
upper_bias <- oc[scenario %in% separated & role == "right" & param == "beta_right" & design %in% c("par", "hyb")]
stopifnot(nrow(upper_bias) == 2 * length(separated))
pct_range <- function(x) sprintf("%.1f–%.1f%%", 100 * min(x), 100 * max(x))
val_range <- function(d) sprintf("%.2f–%.2f", min(d), max(d))
n_reps <- unique(gam$n_reps); stopifnot(length(n_reps) == 1)
gap_vals <- gap[order(arm_f)]

## 4. Assemble and save ---------------------------------------------------------------------
fig <- pA / pB
src <- rbind(
  wd[, .(panel = "a", arm = fit, group = "density", series = "f_model_over_f_kde", x = tau,
         y = ratio_model_over_kde, se = NA_real_)],
  gam[, .(panel = "b", arm = scenario, group = weight_design, series = "gamma_over_truth", x = coverage, y = ratio, se = ratio_se)])
cap <- paste(
  "Figure S7. Recovery of the membership effect under the model and hybrid weight designs.",
  "The model design divides every Stage 2 weight column by the fitted mixture density; the hybrid design divides",
  "the membership (γ) column by the Stage 1 KDE instead (Methods).",
  sprintf("**a**, Fitted mixture density over the Stage 1 KDE at %d quantile levels for testosterone (n = 399,710), unshrunk", n_tau),
  "and at three degrees of shrinkage; each line is one degree of shrinkage and the dashed line marks equality.",
  sprintf("At τ = 0.5 (points), between the modes, the ratios are %s.",
          paste(sprintf("%.2f (%s)", gap_vals$ratio_model_over_kde, c("unshrunk", "0.2", "0.5", "0.7")), collapse = ", ")),
  "**b**, Mean γ estimate over its true value for a variant affecting only membership, in nine simulated scenarios",
  sprintf("(%d replicates each; n = 100,000 unless stated), with 95%% Monte Carlo confidence intervals; the right-hand", n_reps),
  "column gives the coverage of nominal 95% confidence intervals under each design.",
  sprintf("The ratio is %s under the model design and %s under the hybrid design.",
          val_range(gam[weight_design == "model", ratio]), val_range(gam[weight_design == "hybrid", ratio])),
  sprintf("In the %d separated scenarios, the effect of a variant acting only on the narrower, upper component is", length(separated)),
  sprintf("overestimated by %s under the hybrid design and %s under the model design.",
          pct_range(upper_bias[design == "hyb", rel_bias]), pct_range(upper_bias[design == "par", rel_bias])),
  "KDE, kernel density estimate.")
save_figure(fig, "figS07_weight_designs", src, cap, width = 183, height = 180,
            results_used = c("02_mixtures/stage2_mixtures", "02_mixtures/weight_design/sim"))
