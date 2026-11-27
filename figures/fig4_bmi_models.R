#!/usr/bin/env Rscript
# Figure 4: competing BMI models and the observed lead profiles.
# Reproduces the manuscript plotting construction, including the mean-matched
# variance profile in the angle matrix and both reference curves on log BMI.
# Reads aggregate results only. Run from the repository root with Rscript.
suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(data.table)})
.repo <- system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE)
source(file.path(.repo, "figures", "theme.R"))

BASE <- "lognormal"
bl <- as.data.table(read_result("05_bmi/model_profiles", "baseline"))[baseline == BASE]
m0 <- bl$meanlog[1]; s0 <- bl$sdlog[1]; E0 <- exp(m0 + s0^2 / 2)

# ---- A: carrier vs non-carrier (analytic, exaggerated effect for display) --------------------
DELTA <- 2                                   # display mean effect, kg/m^2 (the small-effect profiles are in B)
cm <- log(1 + DELTA / E0)                    # multiplicative: mean + DELTA
s_v <- 1.35 * s0; m_v <- m0 + s0^2 / 2 - s_v^2 / 2   # log-scale variance effect with the mean held fixed
ch <- log(1 + (DELTA / 2) / E0)              # combined: half additive, half multiplicative
mods_a <- c(additive = "Shift", multiplicative = "Proportional\n(shift on log BMI)",
            vqtl = "Pure variance\n(mean-matched)", combined = "Combined\n(shift + proportional)")
x <- seq(14, 52, length.out = 800)
u <- seq(0.02, 0.98, length.out = 200)
dln <- function(x, m, s, shift = 0) ifelse(x > shift, dlnorm(x - shift, m, s), 0)
qln <- function(u, m, s, shift = 0) qlnorm(u, m, s) + shift
dens <- rbind(
  data.table(model = "additive", x, non = dln(x, m0, s0), car = dln(x, m0, s0, DELTA)),
  data.table(model = "multiplicative", x, non = dln(x, m0, s0), car = dln(x, m0 + cm, s0)),
  data.table(model = "vqtl", x, non = dln(x, m0, s0), car = dln(x, m_v, s_v)),
  data.table(model = "combined", x, non = dln(x, m0, s0), car = dln(x, m0 + ch, s0, DELTA / 2)))
qd <- rbind(
  data.table(model = "additive", u, d = qln(u, m0, s0, DELTA) - qln(u, m0, s0)),
  data.table(model = "multiplicative", u, d = qln(u, m0 + cm, s0) - qln(u, m0, s0)),
  data.table(model = "vqtl", u, d = qln(u, m_v, s_v) - qln(u, m0, s0)),
  data.table(model = "combined", u, d = qln(u, m0 + ch, s0, DELTA / 2) - qln(u, m0, s0)))
for (d in list(dens, qd)) d[, model_f := factor(mods_a[model], levels = mods_a)]
dl <- melt(dens, id.vars = c("model_f", "x"), measure.vars = c("non", "car"), variable.name = "who", value.name = "f")
dl[, who := factor(who, levels = c("non", "car"), labels = c("Non-carrier", "Carrier"))]
wcol <- c(`Non-carrier` = "#5b6475", Carrier = "#c4513b")
pA1 <- ggplot(dl, aes(x, f, colour = who, fill = who)) +
  geom_area(position = "identity", alpha = 0.18, linewidth = 0.55, outline.type = "upper") +
  facet_wrap(~model_f, nrow = 1) +
  scale_colour_manual(values = wcol, name = NULL) + scale_fill_manual(values = wcol, name = NULL) +
  scale_x_continuous(breaks = c(20, 30, 40, 50)) +
  labs(x = "BMI (kg/m²)", y = "Density", tag = "a") +
  theme_paper() + legend_inside(0.995, 0.95, c(1, 1)) +
  theme(strip.text = element_text(size = rel(0.9), lineheight = 0.95))
pA2 <- ggplot(qd, aes(u, d)) +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = 2, colour = "grey45") +
  geom_line(linewidth = 0.75, colour = paper_palette[["shift"]], lineend = "round") +
  facet_wrap(~model_f, nrow = 1, scales = "free_y") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", ".25", ".5", ".75", "1"), limits = c(0, 1)) +
  labs(x = "Quantile level tau", y = "Quantile difference\n(kg/m²)") +
  theme_paper() + theme(strip.text = element_blank())
pA <- pA1 / pA2 + plot_layout(heights = c(1, 0.8))

# ---- B: 3 mean-changing curves + 4x4 angle matrix matching Panel a exactly ------------------
mods <- c(additive = "Shift (raw BMI)", log_additive = "Proportional (shift on log BMI)",
          combined = "Combined (shift + proportional)", vqtl_mm = "Pure variance (mean-matched)")
mcol <- c(additive = paper_palette[["shift"]], log_additive = paper_palette[["spread"]],
          combined = paper_palette[["bend"]], vqtl_mm = "#2a9d8f")
lt <- c(additive = "solid", log_additive = "solid", combined = "42", vqtl_mm = "solid")
pr <- as.data.table(read_result("05_bmi/model_profiles", "profiles"))[baseline == BASE & profile %in% c("additive", "log_additive", "combined")]

# Compute exact whitened and quantile-space angles for the 4 Panel-a models (incl. mean-matched log-scale variance)
taus_g <- bl$tau; z_g <- bl$z
s0_bmi <- fread(file.path(.repo, "results", "01_data/lead_profiles", "sigma0_bmi.csv"))[chr == 1]
scols <- grep("^s_", names(s0_bmi), value = TRUE)
S0_bmi <- as.matrix(s0_bmi[, ..scols])
eg <- eigen((S0_bmi + t(S0_bmi)) / 2, symmetric = TRUE)
Sih_bmi <- eg$vectors %*% (pmax(eg$values, 1e-12)^(-0.5) * t(eg$vectors))

X_g <- exp(m0 + s0 * z_g)
# Exact infinitesimal mean-matched log-scale variance profile: d/d(s^2) [exp(m0 + s0^2/2 - s^2/2 + s z)] = X * (z/(2s) - 1/2)
v_mm <- X_g * (z_g / (2 * s0) - 0.5)
P4 <- cbind(
  additive = rep(1, length(taus_g)),
  log_additive = X_g / E0,
  combined = 0.5 * rep(1, length(taus_g)) + 0.5 * (X_g / E0),
  vqtl_mm = v_mm
)
ang_deg <- function(a, b) acos(max(-1, min(1, sum(a * b) / sqrt(sum(a^2) * sum(b^2))))) * 180 / pi
cand <- names(mods)
g <- CJ(a = cand, b = cand)
g[, `:=`(ia = match(a, cand), ib = match(b, cand))]
g <- g[ia < ib]
g[, w := mapply(function(aa, bb) ang_deg(as.numeric(Sih_bmi %*% P4[, aa]), as.numeric(Sih_bmi %*% P4[, bb])), a, b)]
g[, t := mapply(function(aa, bb) ang_deg(P4[, aa], P4[, bb]), a, b)]

pro <- smooth_curve(pr, "tau", "beta_unit_mean", by = "profile")
pro[, model := factor(profile, levels = names(mods))]
pB1 <- ggplot(pro, aes(x, y, colour = model, linetype = model)) +
  geom_hline(yintercept = 1, linewidth = 0.3, colour = "grey60") +
  geom_line(linewidth = 0.75, lineend = "round") +
  scale_colour_manual(values = mcol[1:3], labels = mods[1:3], name = NULL) +
  scale_linetype_manual(values = lt[1:3], labels = mods[1:3], name = NULL) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", ".25", ".5", ".75", "1")) +
  labs(x = "Quantile level tau", y = "Untransformed BMI quantile effect\nper unit mean effect", tag = "b") +
  theme_paper() + legend_inside(0.02, 0.98, c(0, 1)) + theme(legend.key.width = unit(1.4, "lines"))

short <- c(additive = "Shift", log_additive = "Proportional", combined = "Combined", vqtl_mm = "Pure variance")
short2 <- c(additive = "Shift", log_additive = "Proportional", combined = "Combined", vqtl_mm = "Pure\nvariance")
g[, lab := sprintf("%.1f\n(%.1f)", w, t)]
g[, `:=`(a_f = factor(short[a], levels = rev(unname(short))), b_f = factor(short2[b], levels = unname(short2)))]
pB2 <- ggplot(g, aes(b_f, a_f, fill = pmin(w, 90))) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = lab, colour = w > 60), size = 2.1, lineheight = 0.9) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey10"), guide = "none") +
  scale_fill_gradient(low = "#eef3f9", high = "#4f78a8", limits = c(0, 90), guide = "none") +
  scale_x_discrete(limits = unname(short2[-1]), position = "top") + scale_y_discrete(limits = rev(unname(short[-4]))) +
  labs(x = NULL, y = NULL) +
  theme_paper() + theme(panel.grid = element_blank(), panel.grid.major = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(),
                        axis.text.x.top = element_text(lineheight = 0.9))

# ---- C: spaghetti on LOG BMI with BOTH reference curves overlaid & labelled -----------------
lp <- as.data.table(read_result("01_data/lead_profiles", "lead_profiles"))[lead_set == "log_bmi_A_regenie"]
stopifnot(nrow(lp) == 639)
bcols <- grep("^beta_tau_", names(lp), value = TRUE); taus <- tau_grid()
B <- as.matrix(lp[, ..bcols]); S <- B / rowMeans(B)
spg <- data.table(snp_id = rep(lp$snp_id, each = length(taus)), tau = rep(taus, nrow(lp)), scaled = as.vector(t(S)))
med <- spg[, .(median = median(scaled), q25 = quantile(scaled, 0.25), q75 = quantile(scaled, 0.75)), by = tau]
spg_s <- smooth_curve(spg, "tau", "scaled", by = "snp_id", n = 150)
med_s <- merge(smooth_curve(med, "tau", "median"), smooth_curve(med, "tau", "q25")[, .(x, q25 = y)], by = "x")
med_s <- merge(med_s, smooth_curve(med, "tau", "q75")[, .(x, q75 = y)], by = "x")

# Raw-BMI pure shift transformed onto log-BMI scale: (1/q_raw) / mean(1/q_raw)
inv_q <- 1 / bl$q_stage1
raw_shift_on_log <- data.table(tau = bl$tau, y = inv_q / mean(inv_q))
raw_shift_s <- smooth_curve(raw_shift_on_log, "tau", "y")

ref_c <- rbind(
  data.table(x = med_s$x, y = med_s$y, series = "Observed log-BMI median (639 leads)"),
  data.table(x = c(min(med_s$x), max(med_s$x)), y = c(1, 1), series = "Proportional (shift on log BMI)"),
  data.table(x = raw_shift_s$x, y = raw_shift_s$y, series = "Shift on raw BMI")
)
ref_levels <- c("Observed log-BMI median (639 leads)", "Proportional (shift on log BMI)", "Shift on raw BMI")
ref_c[, series := factor(series, levels = ref_levels)]
ref_col <- c("Observed log-BMI median (639 leads)" = paper_palette[["spread"]],
             "Proportional (shift on log BMI)" = "grey20",
             "Shift on raw BMI" = "grey20")
ref_lt  <- c("Observed log-BMI median (639 leads)" = "solid",
             "Proportional (shift on log BMI)" = "dashed",
             "Shift on raw BMI" = "12")
ref_lw  <- c("Observed log-BMI median (639 leads)" = 0.85,
             "Proportional (shift on log BMI)" = 0.50,
             "Shift on raw BMI" = 0.65)

pC <- ggplot(spg_s, aes(x, y, group = snp_id)) +
  geom_line(linewidth = 0.12, alpha = 0.08, colour = paper_palette[["shift"]]) +
  geom_ribbon(data = med_s, aes(x, ymin = q25, ymax = q75), inherit.aes = FALSE, fill = paper_palette[["spread"]], alpha = 0.25) +
  geom_line(data = ref_c, aes(x, y, colour = series, linetype = series, linewidth = series), inherit.aes = FALSE) +
  scale_colour_manual(values = ref_col, name = NULL) +
  scale_linetype_manual(values = ref_lt, name = NULL) +
  scale_linewidth_manual(values = ref_lw, name = NULL) +
  coord_cartesian(ylim = c(-0.3, 2.2)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", ".25", ".5", ".75", "1")) +
  labs(x = "Quantile level tau", y = "Log-BMI quantile effect /\nlead mean effect", tag = "c") +
  theme_paper() + legend_inside(0.02, 0.98, c(0, 1)) + theme(legend.key.width = unit(1.4, "lines"))

# ---- Assemble and export the quantities actually plotted ---------------------------------
fig <- wrap_plots(list(pA, (pB1 | pB2) + plot_layout(widths = c(1.35, 1)), pC),
                  ncol = 1, heights = c(1.15, 1, 0.9))
src <- rbindlist(list(
  qd[, .(panel = "a", series = paste(model, "quantile_difference"), x = u, y = d)],
  dens[, .(panel = "a", series = paste(model, "carrier_density"), x, y = car)],
  dens[model == "additive", .(panel = "a", series = "noncarrier_density", x, y = non)],
  pro[, .(panel = "b", series = as.character(model), x, y)],
  g[, .(panel = "b", series = paste("whitened_angle", a, b), x = NA_real_, y = w)],
  g[, .(panel = "b", series = paste("quantile_angle", a, b), x = NA_real_, y = t)],
  ref_c[, .(panel = "c", series = as.character(series), x, y)],
  med_s[, .(panel = "c", series = "pointwise_interquartile_range", x, y,
             lo = q25, hi = q75)]
), fill = TRUE)
cap <- paste(
  "Figure 4. Competing models of a genetic effect on BMI predict nearly the same quantile profile.",
  "n = 461,320. **a**, Analytic illustration on a log-normal baseline fitted to the Stage 1 BMI quantiles:",
  sprintf("non-carrier and carrier distributions (top) and their quantile difference (bottom) for an allele raising mean BMI by %g kg/m2,", DELTA),
  "exaggerated for display, under an additive shift, a proportional effect, a log-scale variance effect with the",
  "mean held fixed and a combined shift and proportional effect. **b**, Shift, proportional and combined quantile-effect",
  "profiles per unit mean effect (left), and pairwise angles for all four panel-a models in the whitened metric,",
  "with quantile-space angles in brackets (right). **c**, Stage 1 log-BMI profiles of 639 lead variants, each divided",
  "by its own mean effect (blue); the red line and band show the pointwise median and interquartile range.",
  "Reference curves show a shift on log BMI and a shift on raw BMI transported to log BMI.",
  "The simulation-calibrated, equal-rank comparison rejects a raw-scale additive shift at 570 of 639 leads (89.2%)",
  "on log BMI and the reverse comparison at one lead (Methods).")
save_figure(fig, "fig4_bmi_models", src, cap, width = 183, height = 170,
            results_used = c("05_bmi/model_profiles/baseline", "05_bmi/model_profiles/profiles",
                             "01_data/lead_profiles/lead_profiles", "01_data/lead_profiles/sigma0_bmi",
                             "05_bmi/model_separation/lead_pair_power"))
