#!/usr/bin/env Rscript
# Supplementary Figure S10: moment and principal-component descriptions of BMI lead profiles
#
# What: two ways to describe how a variant moves the BMI distribution. Moments use shapes fixed
#   in advance (shift, spread, skew as L-moment weight curves); principal components use the
#   shapes the 639 lead profiles themselves choose. Each lead profile is divided by the lead's
#   mean effect before comparison.
#   a: the L-moment weight curves.
#   b: each lead's effect at tau 0.89 against tau 0.11, with moment axes and principal axes.
#   c: whitened principal-component spectrum against a parallel-analysis null.
#   d: what is left of every lead's profile after fitting PC1, PC1 + PC2, shift, shift + spread.
# Paper: Supplementary Fig. S10
# Input: results/01_data/lead_profiles/lead_profiles.csv,
#   results/06_traits/moments/{pca_loadings,pca_spectrum,pca_parallel}.csv
# Output: paper_assets/figS10_moments_pcs.{png,pdf}, _source.csv, _caption.md
# Runs: locally from results/ in about 5 seconds

suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(data.table)})
source(file.path(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE), "figures", "theme.R"))

## 1. Lead profiles scaled by their mean effect -------------------------------------------
lp <- as.data.table(read_result("01_data/lead_profiles", "lead_profiles"))[lead_set == "log_bmi_A_regenie"]
stopifnot(nrow(lp) == 639)
beta_cols <- grep("^beta_tau_", names(lp), value = TRUE); tau <- tau_grid(); stopifnot(length(beta_cols) == length(tau))
z <- qnorm(tau)
S <- as.matrix(lp[, ..beta_cols]) / lp$delta_mu      # effect / mean effect
U <- S / sqrt(rowSums(S^2))                           # unit-length profiles
# PC loadings mapped to tau space (unwhitened), signed so PC1 is positive and PC2 rises.
loadings <- as.data.table(read_result("06_traits/moments", "pca_loadings"))[arm == "log_bmi"]
unit_len <- function(x) x / sqrt(sum(x^2))
v1 <- unit_len(loadings[pc == 1, loading_tau]); v2 <- unit_len(loadings[pc == 2, loading_tau])
if (sum(v1) < 0) v1 <- -v1; if (tail(v2, 1) < 0) v2 <- -v2
shift_curve <- unit_len(rep(1, length(tau))); spread_curve <- unit_len(z)   # Cornish-Fisher shift and spread
proj <- function(Y, X) { X <- as.matrix(X); Y %*% X %*% solve(crossprod(X), t(X)) }   # least-squares fit of rows of Y
BL <- paper_palette[["shift"]]; CO <- paper_palette[["spread"]]; GR <- "grey60"; TE <- "#1D9E75"
tau_axis <- scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", ".5", "1"))

## 2. Panel a: L-moment weight curves -----------------------------------------------------
wl <- c("Shift", "Spread", "Skew")
weights <- data.table(tau = rep(tau, 3), w = c(rep(1, length(tau)), 2 * tau - 1, 6 * tau^2 - 6 * tau + 1),
                      what = factor(rep(wl, each = length(tau)), wl))
pA <- ggplot(smooth_curve(weights, "tau", "w", by = "what"), aes(x, y, colour = what)) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_line(linewidth = 0.9, lineend = "round") +
  scale_colour_manual(values = setNames(c(BL, CO, TE), wl), name = NULL) + tau_axis +
  labs(x = expression("Quantile level " * tau), y = "L-moment weight", tag = "a") +
  theme_paper() + legend_inside(0.62, 0.02, c(0.5, 0))

## 3. Panel b: two percentiles ------------------------------------------------------------
i1 <- which.min(abs(tau - 0.107)); i2 <- which.min(abs(tau - 0.893))
P <- data.table(x = S[, i1], y = S[, i2])
e <- eigen(crossprod(as.matrix(P)))$vectors[, 1]; if (e[1] < 0) e <- -e   # principal axis through the origin
mu <- colMeans(P); ang <- atan2(e[2], e[1]) * 180 / pi - 45               # its angle above the shift axis
pB <- ggplot(P, aes(x, y)) + geom_point(colour = GR, alpha = 0.45, size = 0.6, stroke = 0) +
  geom_abline(slope = 1, intercept = 0, colour = BL, linetype = 2, linewidth = 0.6) +    # pure shift
  geom_abline(slope = -1, intercept = 2, colour = BL, linetype = 2, linewidth = 0.6) +   # pure spread
  geom_abline(slope = e[2] / e[1], intercept = 0, linewidth = 0.7) +
  geom_abline(slope = -e[1] / e[2], intercept = mu[2] + e[1] / e[2] * mu[1], linewidth = 0.7) +
  coord_equal(xlim = c(-0.2, 2.5), ylim = c(-0.35, 2.6)) +
  labs(x = bquote("Effect at " * tau == .(sprintf("%.2f", tau[i1])) * " / mean effect"),
       y = bquote("Effect at " * tau == .(sprintf("%.2f", tau[i2])) * " / mean effect"), tag = "b") + theme_paper()

## 4. Panel c: principal-component spectrum -----------------------------------------------
spectrum <- as.data.table(read_result("06_traits/moments", "pca_spectrum"))[arm == "log_bmi" & pc <= 8]
parallel <- as.data.table(read_result("06_traits/moments", "pca_parallel"))[arm == "log_bmi"]
null95 <- parallel[, .(nul = null_p95_pct[1]), by = .(pc = component_tested)]
n_comp <- parallel[, n_components[1]]
pC <- ggplot(spectrum, aes(pc, pct)) + geom_line(colour = "grey55", linewidth = 0.4) +
  geom_point(data = null95[pc <= 8], aes(pc, nul), shape = 95, size = 6, colour = "grey60") +
  geom_point(aes(colour = pc <= n_comp), size = 1.8) +
  scale_colour_manual(values = c(`TRUE` = CO, `FALSE` = "grey35"), guide = "none") +
  scale_y_log10() + scale_x_continuous(breaks = 1:8) +
  labs(x = "Principal component", y = "Shape variance explained (%)", tag = "c") + theme_paper()

## 5. Panel d: leftovers after each fit ---------------------------------------------------
# Unwhitened unit-length profiles, so the fractions are not comparable with panel c.
fits <- list(`PC1` = list(v1, "black"), `PC1 + PC2` = list(cbind(v1, v2), "black"),
             `shift` = list(shift_curve, BL), `shift + spread` = list(cbind(shift_curve, spread_curve), BL))
fit_lab <- c(`PC1` = "PC1", `PC1 + PC2` = "PC1 + PC2", shift = "Shift", `shift + spread` = "Shift + spread")
left <- sapply(fits, function(f) sum((U - proj(U, f[[1]]))^2) / sum(U^2))
leftover <- rbindlist(lapply(names(fits), function(n) {
  R <- U - proj(U, fits[[n]][[1]])
  data.table(fit = n, id = rep(seq_len(nrow(R)), each = length(tau)), tau = rep(tau, nrow(R)), v = as.vector(t(R)))
}))
leftover[, fit_f := factor(fit_lab[fit], levels = fit_lab)]
pD <- ggplot(leftover, aes(tau, v, group = id, colour = fit_f)) + geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.3) +
  geom_line(alpha = 0.06, linewidth = 0.25) + coord_cartesian(ylim = c(-0.25, 0.3)) + tau_axis +
  facet_wrap(~fit_f, nrow = 1) +
  scale_colour_manual(values = setNames(sapply(fits, `[[`, 2), fit_lab), guide = "none") +
  labs(x = expression("Quantile level " * tau), y = "Leftover", tag = "d") + theme_paper()

## 6. Assemble and save ---------------------------------------------------------------------
fig <- wrap_plots(list((pA | pB | pC) + plot_layout(widths = c(1, 1, 1)), pD), ncol = 1, heights = c(1.3, 1))
src <- rbind(
  P[, .(panel = "b", series = lp$snp_id, x, y)],
  spectrum[, .(panel = "c", series = "observed_pct", x = pc, y = pct)], null95[, .(panel = "c", series = "null_p95_pct", x = pc, y = nul)],
  data.table(panel = "d", series = paste(names(fits), "fraction_left"), x = NA_real_, y = unname(left)))
cap <- paste(
  "Figure S10. Moment and principal-component descriptions of BMI lead-variant profiles.",
  sprintf("Log BMI; %d lead variants (BMI-increasing allele), each Stage 1 profile divided by the lead's mean effect;", nrow(lp)),
  "n = 461,320.",
  "**a**, L-moment weight curves for shift (1), spread (2τ − 1) and skew (6τ² − 6τ + 1): the average of a profile",
  "multiplied by a curve gives that moment.",
  sprintf("**b**, Effect at τ = %.2f against effect at τ = %.2f, each over the mean effect; each point is one lead.", tau[i2], tau[i1]),
  "Blue dashed lines are the moment axes (pure shift through (1, 1) and pure spread orthogonal to it) and black lines",
  sprintf("the principal axes of the cloud; the first lies %.0f° above the shift axis.", ang),
  "**c**, Share of shape variance explained by each principal component of the whitened unit-length profiles (log scale)",
  sprintf("against the 95th percentile of a parallel-analysis null (ticks); red, the %d components above the null.", n_comp),
  "**d**, What is left of each lead's unit-length profile after fitting PC1, PC1 + PC2, shift, and shift + spread",
  sprintf("(%s of total squared length; unwhitened, so not comparable with the whitened shares).",
          paste(sprintf("%.1f%%", 100 * left), collapse = ", ")),
  "The second component is compared across traits in Fig. 5b.")
save_figure(fig, "figS10_moments_pcs", src, cap, width = 183, height = 125,
            results_used = c("01_data/lead_profiles", "06_traits/moments"))
