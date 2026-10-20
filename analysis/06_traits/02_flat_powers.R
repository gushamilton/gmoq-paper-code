#!/usr/bin/env Rscript
# "Flat" Box-Cox powers: the scale on which genetic effects would be pure shifts
#
# What: For each trait, every lead profile is moved to the raw measurement scale and then to
#   each Box-Cox power lambda in -5 ... 5 by exact change of scale. A lead's flat power is the
#   lambda at which its converted profile has zero L-moment spread. The pooled flat power is the
#   Box-Cox curve that best matches the trait's implied natural scale g(u) = int dQ / pbar, where
#   pbar is the precision-weighted mean lead profile on the raw scale: on g the average lead is a
#   pure shift. Supplementary Table S5 sets these against the normalising power of the trait.
# Paper: Supplementary Table S5; Methods "Change of scale"; Discussion (flat power of BMI leads)
# Input: aligned lead-variant profiles and per-chromosome pooled covariance
#   (results/01_data/lead_profiles), baseline quantiles (results/01_data/baseline_quantiles)
# Output: results/06_traits/scale/{flat_lambda_summary,flat_lambda_leads,implied_scale,
#   implied_scale_comparison}.csv
# Runs: locally from results/ in a few seconds

repo <- Sys.getenv("GMOQ_REPO", system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE))
source(file.path(repo, "R", "load.R"))

cfg <- read_config()
u <- tau_grid(cfg); n_tau <- length(u); z <- qnorm(u)
cw <- trap_weights(u)                     # trapezoid weights: delta_mu = sum cw * beta
w2 <- cw * (2 * u - 1)                    # second L-moment (spread) weights
lambda_grid <- round(seq(-500, 500) / 100, 2)

lead_profiles <- as.data.table(read_result("01_data/lead_profiles", "lead_profiles"))
baseline <- as.data.table(read_result("01_data/baseline_quantiles", "baseline_quantiles"))
baseline_quantiles <- function(trait_name) {
  x <- baseline[baseline$trait == trait_name]
  as.numeric(x$q[match(round(u, 3), round(as.numeric(x$tau), 3))])
}

# Traits: raw BMI, weight and the survey traits, one row each (log BMI is BMI on another scale).
survey <- cfg$cohorts$survey$traits
survey_traits <- names(survey)[vapply(survey, function(x) !identical(x$fig6, FALSE) &&
                                        !is.null(x[["n", exact = TRUE]]), logical(1))]
traits <- data.table(trait = c("bmi", "weight", survey_traits),
                     lead_set = c("bmi_A_regenie", "weight_A_regenie", paste0(survey_traits, "_A_deltamu")))

## 1. Lead profiles on the raw measurement scale
# A log-analysed trait is converted back exactly: beta_raw(u) = Q_raw(u) beta_log(u) and
# Sigma_raw = D Sigma0 D with D = diag(Q_raw(u)).
raw_profiles <- function(trait_name, lead_set) {
  leads <- lead_set_matrix(lead_profiles, lead_set, u)
  q_analysed <- baseline_quantiles(trait_name)
  is_log <- grepl("^log_", trait_name)
  q_raw <- if (is_log) exp(q_analysed) else q_analysed
  B_raw <- if (is_log) sweep(leads$B, 2, q_raw, `*`) else leads$B
  sigma0 <- p05_sigma0(leads$source, u)
  sigma0_raw <- if (is_log) lapply(sigma0, function(S) change_scale(matrix(0, 1, n_tau), S, q_raw)$Sigma0) else sigma0
  list(meta = leads$meta, q_raw = q_raw, B_raw = B_raw, sigma0_raw = sigma0_raw)
}

# Root nearest zero of the spread T(lambda) along the lambda grid, by linear interpolation;
# -Inf / +Inf when the spread keeps one sign over the whole grid.
root_nearest_zero <- function(t) {
  k <- which(diff(sign(t)) != 0)
  if (length(k)) {
    k <- k[which.min(abs(lambda_grid[k]))]
    lambda_grid[k] - t[k] * (lambda_grid[k + 1] - lambda_grid[k]) / (t[k + 1] - t[k])
  } else if (t[1] > 0) -Inf else Inf
}

implied <- list(); implied_comparison <- list(); flat_leads <- list(); flat_summary <- list()
for (i in seq_len(nrow(traits))) {
  tr <- traits$trait[i]
  rp <- raw_profiles(tr, traits$lead_set[i])
  B <- rp$B_raw; q_raw <- rp$q_raw; d <- rp$meta$d_resid; chr <- as.character(rp$meta$chr)
  positive <- all(q_raw > 0)              # Box-Cox needs a positive trait

  ## 2. Implied natural scale and the pooled flat power
  # Align each lead to its effect-increasing allele, average with weights d (precision), and
  # integrate dQ / pbar (trapezoid): on the scale g(Q) every quantile moves by the same amount.
  sgn <- sign(as.numeric(B %*% cw)); sgn[sgn == 0] <- 1
  pbar <- as.numeric(colSums(d * sgn * B) / sum(d))
  g <- cumsum(c(0, diff(q_raw) / ((pbar[-1] + pbar[-n_tau]) / 2)))
  fit_r <- vapply(lambda_grid, function(l) if (!positive) NA_real_ else suppressWarnings(cor(boxcox(q_raw, l), g)), 0)
  lambda_pooled <- if (all(is.na(fit_r))) NA_real_ else lambda_grid[which.max(fit_r)]
  implied[[tr]] <- data.table(trait = tr, tau = u, q_raw = q_raw, pbar_raw = pbar, g = g, z = z)
  implied_comparison[[tr]] <- data.table(
    trait = tr, comparator = c("rank_inverse_normal", "log", "identity", "boxcox_pooled_flat"),
    pearson_r = c(cor(g, z), if (positive) cor(g, log(q_raw)) else NA, cor(g, q_raw),
                  if (is.na(lambda_pooled)) NA else max(fit_r, na.rm = TRUE)))

  ## 3. Per-lead flat power
  # On the Box-Cox scale lambda the profile is Q_raw^(lambda - 1) beta_raw, so its L-moment
  # spread is T(lambda) = sum w2 Q_raw^(lambda - 1) beta_raw (sign-aligned). The flat power is
  # the root of T nearest zero. "Strong" leads have |raw-scale mean-effect z| > 10.
  mean_z <- sgn * as.numeric(B %*% cw) /
    sqrt(vapply(seq_along(d), function(k) functional_cov(cw, cw, rp$sigma0_raw[[chr[k]]], d[k]), 0))
  spread_by_lambda <- if (!positive) NULL else
    vapply(lambda_grid, function(l) as.numeric(B %*% (w2 * q_raw^(l - 1))) * sgn, numeric(nrow(B)))
  lambda_flat <- if (is.null(spread_by_lambda)) rep(NA_real_, nrow(B)) else apply(spread_by_lambda, 1, root_nearest_zero)
  clipped <- pmin(pmax(lambda_flat, -5), 5)  # no root: counted at the edge of the search
  strong <- abs(mean_z) > 10
  q_or_na <- function(x, p) if (length(x)) quantile(x, p, na.rm = TRUE, names = FALSE) else NA_real_
  flat_leads[[tr]] <- data.table(trait = tr, snp_id = rp$meta$snp_id, chr = rp$meta$chr,
                                 raw_mean_z = mean_z, lambda_flat_v = lambda_flat)
  flat_summary[[tr]] <- data.table(
    trait = tr, n_leads = length(lambda_flat), n_strong = sum(strong),
    lambda_v_median = median(clipped, na.rm = TRUE), lambda_v_q25 = q_or_na(clipped, 0.25),
    lambda_v_q75 = q_or_na(clipped, 0.75),
    lambda_v_median_strong = if (any(strong)) median(clipped[strong], na.rm = TRUE) else NA_real_,
    lambda_v_q25_strong = if (any(strong)) q_or_na(clipped[strong], 0.25) else NA_real_,
    lambda_v_q75_strong = if (any(strong)) q_or_na(clipped[strong], 0.75) else NA_real_,
    frac_no_root = mean(!is.finite(lambda_flat)))
}

## 4. Write
out <- "06_traits/scale"
write_result(rbindlist(flat_summary), out, "flat_lambda_summary")
write_result(rbindlist(flat_leads), out, "flat_lambda_leads")
write_result(rbindlist(implied), out, "implied_scale")
write_result(rbindlist(implied_comparison), out, "implied_scale_comparison")
