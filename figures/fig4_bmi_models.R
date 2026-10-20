#!/usr/bin/env Rscript
# Figure 4: competing models of a genetic effect on BMI predict nearly the same quantile profile
#
# What: an additive effect, a multiplicative (log-additive) effect, a variance effect and
#   their combinations are different mechanisms, but on a BMI-like right-skewed baseline
#   they predict almost the same quantile-effect profile. The figure shows this analytically
#   and then shows that the 639 BMI lead variants share one common tilted profile.
#   a: carrier and non-carrier distributions and their quantile difference under four models.
#   b: the model profiles per unit mean effect, and the pairwise angles between them.
#   c: the Stage 1 profiles of the 639 BMI leads, each divided by its own mean effect.
# Paper: Fig. 4
# Input: results/05_bmi/model_profiles/{baseline,profiles,angles_pairwise}.csv,
#   results/01_data/lead_profiles/lead_profiles.csv
# Output: paper_assets/fig4_bmi_models.{png,pdf}, _source.csv, _caption.md
# Runs: locally from results/ in about 5 seconds

suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(data.table)})
source(file.path(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE), "figures", "theme.R"))

tau_axis <- scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", ".25", ".5", ".75", "1"))

## 1. Load ----------------------------------------------------------------------------------
# The log-normal baseline fitted to the Stage 1 BMI quantiles.
baseline <- as.data.table(read_result("05_bmi/model_profiles", "baseline"))[baseline == "lognormal"]
profiles <- as.data.table(read_result("05_bmi/model_profiles", "profiles"))[baseline == "lognormal"]
angles <- as.data.table(read_result("05_bmi/model_profiles", "angles_pairwise"))[baseline == "lognormal"]
lead_prof <- as.data.table(read_result("01_data/lead_profiles", "lead_profiles"))[lead_set == "log_bmi_A_regenie"]
stopifnot(nrow(lead_prof) == 639)

## 2. Panel a: carrier against non-carrier, analytically -------------------------------------
m0 <- baseline$meanlog[1]; s0 <- baseline$sdlog[1]; E0 <- exp(m0 + s0^2 / 2)
DELTA <- 2                                   # display mean effect (kg/m^2), exaggerated for display
cm <- log(1 + DELTA / E0)                    # multiplicative effect with the same mean change
s_v <- 1.35 * s0; m_v <- m0 + s0^2 / 2 - s_v^2 / 2   # log-scale variance effect, mean held fixed
ch <- log(1 + (DELTA / 2) / E0)              # combined: half additive, half multiplicative
mods_a <- c(additive = "Additive", multiplicative = "Multiplicative (log-additive)",
            vqtl = "vQTL (log-scale variance,\nmean matched)", combined = "Combined\n(additive + multiplicative)")
x <- seq(14, 52, length.out = 800)
u <- seq(0.02, 0.98, length.out = 200)
dln <- function(x, m, s, shift = 0) ifelse(x > shift, dlnorm(x - shift, m, s), 0)   # shifted log-normal density
qln <- function(u, m, s, shift = 0) qlnorm(u, m, s) + shift                         # and quantile function
dens <- rbind(
  data.table(model = "additive", x, non = dln(x, m0, s0), car = dln(x, m0, s0, DELTA)),
  data.table(model = "multiplicative", x, non = dln(x, m0, s0), car = dln(x, m0 + cm, s0)),
  data.table(model = "vqtl", x, non = dln(x, m0, s0), car = dln(x, m_v, s_v)),
  data.table(model = "combined", x, non = dln(x, m0, s0), car = dln(x, m0 + ch, s0, DELTA / 2)))
qdiff <- rbind(
  data.table(model = "additive", u, d = qln(u, m0, s0, DELTA) - qln(u, m0, s0)),
  data.table(model = "multiplicative", u, d = qln(u, m0 + cm, s0) - qln(u, m0, s0)),
  data.table(model = "vqtl", u, d = qln(u, m_v, s_v) - qln(u, m0, s0)),
  data.table(model = "combined", u, d = qln(u, m0 + ch, s0, DELTA / 2) - qln(u, m0, s0)))
for (d in list(dens, qdiff)) d[, model_f := factor(mods_a[model], levels = mods_a)]
dens_long <- melt(dens, id.vars = c("model_f", "x"), measure.vars = c("non", "car"), variable.name = "who", value.name = "f")
dens_long[, who := factor(who, levels = c("non", "car"), labels = c("Non-carrier", "Carrier"))]
who_col <- c(`Non-carrier` = "#5b6475", Carrier = "#c4513b")
pA_dens <- ggplot(dens_long, aes(x, f, colour = who, fill = who)) +
  geom_area(position = "identity", alpha = 0.18, linewidth = 0.55, outline.type = "upper") +
  facet_wrap(~model_f, nrow = 1) +
  scale_colour_manual(values = who_col, name = NULL) + scale_fill_manual(values = who_col, name = NULL) +
  scale_x_continuous(breaks = c(20, 30, 40, 50)) +
  labs(x = "BMI (kg/m²)", y = "Density", tag = "a") +
  theme_paper() + legend_inside(0.995, 0.95, c(1, 1)) +
  theme(strip.text = element_text(size = rel(0.9), lineheight = 0.95))
pA_qdiff <- ggplot(qdiff, aes(u, d)) +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = 2, colour = "grey45") +
  geom_line(linewidth = 0.75, colour = paper_palette[["shift"]], lineend = "round") +
  facet_wrap(~model_f, nrow = 1, scales = "free_y") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", ".25", ".5", ".75", "1"), limits = c(0, 1)) +
  labs(x = "Quantile level tau", y = "Quantile difference\n(kg/m^2)") +
  theme_paper() + theme(strip.text = element_blank())
pA <- pA_dens / pA_qdiff + plot_layout(heights = c(1, 0.8))

## 3. Panel b: model profiles and the angles between them ----------------------------------
mods <- c(additive = "Additive", log_additive = "Log-additive (multiplicative)", combined = "Combined",
          vqtl_raw = "Raw-scale vQTL (constant CV)", vqtl_log = "Log-scale vQTL")
mod_col <- c(additive = paper_palette[["shift"]], log_additive = paper_palette[["spread"]], combined = paper_palette[["bend"]],
             vqtl_raw = "#e4a02e", vqtl_log = "#2a9d8f")
mod_lty <- c(additive = "solid", log_additive = "solid", combined = "42", vqtl_raw = "13", vqtl_log = "solid")
profiles <- profiles[profile %in% names(mods)]
# The log-scale vQTL profile is left out of the overlay (it is in the angle matrix).
prof_curves <- smooth_curve(profiles[profile != "vqtl_log"], "tau", "beta_unit_mean", by = "profile")
prof_curves[, model := factor(profile, levels = names(mods))]
pB_prof <- ggplot(prof_curves, aes(x, y, colour = model, linetype = model)) +
  geom_hline(yintercept = 1, linewidth = 0.3, colour = "grey60") +
  geom_line(linewidth = 0.75, lineend = "round") +
  scale_colour_manual(values = mod_col, labels = mods, name = NULL) +
  scale_linetype_manual(values = mod_lty, labels = mods, name = NULL) +
  tau_axis +
  labs(x = "Quantile level tau", y = "Quantile effect per unit mean effect", tag = "b") +
  theme_paper() + legend_inside(0.02, 0.98, c(0, 1)) + theme(legend.key.width = unit(1.4, "lines"))

# Upper triangle of pairwise angles: whitened metric, with the angle in tau space in brackets.
cand <- names(mods)
pair_key <- function(a, b) paste(pmin(a, b), pmax(a, b))
g <- CJ(a = cand, b = cand); g[, k := pair_key(a, b)]
angles[, k := pair_key(profile_a, profile_b)]
g <- merge(g, angles[metric == "whitened", .(k, w = angle_deg)], by = "k", all.x = TRUE)
g <- merge(g, angles[metric == "tau", .(k, t = angle_deg)], by = "k", all.x = TRUE)
short <- c(additive = "Additive", log_additive = "Log-additive", combined = "Combined", vqtl_raw = "vQTL raw", vqtl_log = "vQTL log")
short2 <- c(additive = "Additive", log_additive = "Log-\nadditive", combined = "Combined", vqtl_raw = "vQTL\nraw", vqtl_log = "vQTL\nlog")
g[, `:=`(ia = match(a, cand), ib = match(b, cand))]
g <- g[ia < ib]
g[, lab := sprintf("%.1f\n(%.1f)", w, t)]
g[, `:=`(a_f = factor(short[a], levels = rev(unname(short))), b_f = factor(short2[b], levels = unname(short2)))]
pB_angles <- ggplot(g, aes(b_f, a_f, fill = pmin(w, 90))) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = lab, colour = w > 60), size = 2.1, lineheight = 0.9) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey10"), guide = "none") +
  scale_fill_gradient(low = "#eef3f9", high = "#4f78a8", limits = c(0, 90), guide = "none") +
  scale_x_discrete(limits = unname(short2[-1]), position = "top") + scale_y_discrete(limits = rev(unname(short[-5]))) +
  labs(x = NULL, y = NULL) +
  theme_paper() + theme(panel.grid = element_blank(), panel.grid.major = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(),
                        axis.text.x.top = element_text(lineheight = 0.9))

## 4. Panel c: the 639 BMI lead profiles, each scaled by its mean effect ---------------------
beta_cols <- grep("^beta_tau_", names(lead_prof), value = TRUE); taus <- tau_grid()
B <- as.matrix(lead_prof[, ..beta_cols]); scaled <- B / rowMeans(B)
spag <- data.table(snp_id = rep(lead_prof$snp_id, each = length(taus)), tau = rep(taus, nrow(lead_prof)), scaled = as.vector(t(scaled)))
med <- spag[, .(median = median(scaled), q25 = quantile(scaled, 0.25), q75 = quantile(scaled, 0.75)), by = tau]
spag_s <- smooth_curve(spag, "tau", "scaled", by = "snp_id", n = 150)
med_s <- merge(smooth_curve(med, "tau", "median"), smooth_curve(med, "tau", "q25")[, .(x, q25 = y)], by = "x")
med_s <- merge(med_s, smooth_curve(med, "tau", "q75")[, .(x, q75 = y)], by = "x")
pC <- ggplot(spag_s, aes(x, y, group = snp_id)) +
  geom_line(linewidth = 0.12, alpha = 0.09, colour = paper_palette[["shift"]]) +
  geom_hline(yintercept = 1, linewidth = 0.35, linetype = 2, colour = "grey25") +
  geom_ribbon(data = med_s, aes(x, ymin = q25, ymax = q75), inherit.aes = FALSE, fill = paper_palette[["spread"]], alpha = 0.25) +
  geom_line(data = med_s, aes(x, y), inherit.aes = FALSE, colour = paper_palette[["spread"]], linewidth = 0.85) +
  coord_cartesian(ylim = c(-0.5, 2.4)) +
  tau_axis +
  labs(x = "Quantile level tau", y = "Effect / lead's mean effect", tag = "c") +
  theme_paper()

## 5. Assemble and save ---------------------------------------------------------------------
fig <- wrap_plots(list(pA, (pB_prof | pB_angles) + plot_layout(widths = c(1.35, 1)), pC),
                  ncol = 1, heights = c(1.15, 1, 0.9))
src <- rbind(
  qdiff[, .(panel = "A", series = paste(model, "quantile_difference"), x = u, y = d)],
  dens[seq(1, .N, by = 4), .(panel = "A", series = paste(model, "carrier_density"), x, y = car)],
  dens[model == "additive"][seq(1, .N, by = 4), .(panel = "A", series = "noncarrier_density", x, y = non)],
  profiles[, .(panel = "B", series = profile, x = tau, y = beta_unit_mean)],
  angles[, .(panel = "B", series = paste(metric, profile_a, profile_b), x = NA_real_, y = angle_deg)],
  med[, .(panel = "C", series = "median_scaled_profile", x = tau, y = median)])
cap <- paste(
  "Figure 4. Competing models of a genetic effect on BMI predict nearly the same quantile profile.",
  "n = 461,320. **a**, Analytic illustration on a log-normal baseline fitted to the Stage 1 BMI quantiles:",
  sprintf("non-carrier and carrier distributions (top) and their quantile difference (bottom) for an allele raising mean BMI by %g kg/m2,", DELTA),
  "exaggerated for display, under an additive shift, a multiplicative effect, a log-scale variance effect with the",
  "mean held fixed and a combined additive and multiplicative effect. **b**, The same models as quantile-effect",
  "profiles per unit mean effect (left) and the pairwise angles between them in the whitened metric, with the angle in",
  "quantile space in brackets (right). A raw-scale variance effect at constant coefficient of variation is exactly the",
  "multiplicative profile on this baseline. **c**, Stage 1 profiles of the 639 BMI lead variants, each divided by its",
  "own mean effect (blue); the red line and band are the pointwise median and interquartile range and the dashed line",
  "is a pure shift. Models are fitted to the same Stage 1 output by ordinary least squares with the pooled Stage 1",
  "covariance; a purely additive profile is rejected in favour of a multiplicative one at 89% of leads (90% on the raw scale), and the reverse",
  "at one lead (Methods).")
save_figure(fig, "fig4_bmi_models", src, cap, width = 183, height = 170,
            results_used = c("05_bmi/model_profiles", "01_data/lead_profiles"))
