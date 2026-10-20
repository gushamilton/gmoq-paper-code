#!/usr/bin/env Rscript
# Quantile-effect profiles predicted by competing models of a BMI variant
#
# What: on a log-normal baseline fitted to the 45 Stage 1 BMI quantiles, works out the
#   quantile profile that each candidate model predicts for a variant (additive,
#   multiplicative, their combination, a raw-scale and a log-scale variance effect), each
#   scaled to raise mean BMI by one unit. It then measures how far apart these profiles are:
#   the angle between every pair, in the whitened metric (the pooled Stage 1 covariance
#   removed, so angles reflect what the data can distinguish) and in plain quantile space,
#   and the size of effect a variant needs before a test can tell one model from another.
# Paper: Fig. 4a,b
# Input: Stage 1 BMI quantiles and summary (results/01_data/baseline_quantiles) and the
#   pooled Stage 1 covariance of raw BMI per chromosome (results/01_data/lead_profiles)
# Output: results/05_bmi/model_profiles/{baseline,weight_gates,profiles,angles_pairwise,
#   angles_to_span,separation,sigma0_check}.csv
# Runs: locally from results/ in a few seconds

repo <- Sys.getenv("GMOQ_REPO", unset = NA)
if (is.na(repo)) repo <- system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE)
source(file.path(repo, "R", "load.R"))
source(file.path(repo, "analysis", "05_bmi", "bmi_common.R"))
cfg <- read_config()
u <- tau_grid(cfg); z <- qnorm(u); n_tau <- length(u)
OUT <- "05_bmi/model_profiles"

## 1. Inputs: Stage 1 quantiles of BMI and the pooled covariance ---------------------------
baseline_q <- read_result("01_data/baseline_quantiles", "baseline_quantiles")
baseline_summary <- read_result("01_data/baseline_quantiles", "baseline_summary")
bmi_summary <- baseline_summary[baseline_summary$trait == "bmi" & baseline_summary$group == "main", ]
q_stage1 <- stage1_quantiles(baseline_q, "bmi", u)$q

# One pooled Sigma0 per chromosome; the metric used here is their average, which differs
# from every chromosome's own matrix by less than 0.1% (recorded in sigma0_check).
sigma0_by_chr <- read_lead_sigma0("bmi", u)[as.character(1:22)]
Sigma0 <- Reduce(`+`, sigma0_by_chr) / length(sigma0_by_chr)
sigma0_check <- data.table(chr = 1:22,
  max_rel_diff_vs_mean = vapply(sigma0_by_chr, function(S) max(abs(S - Sigma0)) / max(abs(Sigma0)), 0),
  diag_min = vapply(sigma0_by_chr, function(S) min(diag(S)), 0),
  diag_max = vapply(sigma0_by_chr, function(S) max(diag(S)), 0))
write_result(sigma0_check, OUT, "sigma0_check")

## 2. Log-normal baseline ---------------------------------------------------------------
# log q(u) = m + s z_u, so a least-squares line of log Q(u) on z_u gives m and s (c = 0).
line_fit <- lm(log(q_stage1) ~ z)
b <- c(c = 0, m = unname(coef(line_fit)[1]), s = unname(coef(line_fit)[2]))
X <- exp(b[["m"]] + b[["s"]] * z)                         # q(u) - c
E <- exp(b[["m"]] + b[["s"]]^2 / 2)                       # mean of the log-normal part

## 3. Weights and their check --------------------------------------------------------------
# W[u, p] = dq_u / dtheta_p, the parametric-model-density weight -dF(q; theta)/dtheta / f(q; theta):
#   additive (c): 1;  log-additive (m): X;  log-scale variance (s^2): X z / (2 s).
W <- cbind(additive = 1, log_additive = X, sigma2 = X * z / (2 * b[["s"]]))
# The same weights in density-ratio form, by finite differences of the log-normal CDF. The
# two must agree; this confirms the analytic derivatives.
q_model <- b[["c"]] + X
cdf <- function(shift, meanlog, var_log) plnorm(q_model - shift, meanlog, sqrt(var_log))
f_model <- dlnorm(q_model - b[["c"]], b[["m"]], b[["s"]])
h <- 1e-6; s2 <- b[["s"]]^2
W_density_ratio <- cbind(
  additive = -(cdf(b[["c"]] + h, b[["m"]], s2) - cdf(b[["c"]] - h, b[["m"]], s2)) / (2 * h) / f_model,
  log_additive = -(cdf(b[["c"]], b[["m"]] + h, s2) - cdf(b[["c"]], b[["m"]] - h, s2)) / (2 * h) / f_model,
  sigma2 = -(cdf(b[["c"]], b[["m"]], s2 + h) - cdf(b[["c"]], b[["m"]], s2 - h)) / (2 * h) / f_model)
rel_diff <- apply(abs(W - W_density_ratio), 2, max) / apply(abs(W), 2, max)
weight_gates <- data.table(baseline = "lognormal", param = colnames(W),
                           max_rel_diff_analytic_vs_density_ratio = rel_diff, pass = rel_diff < 1e-5)
stopifnot(all(weight_gates$pass))

baseline_table <- data.table(baseline = "lognormal", tau = u, z = z, q_model = q_model, q_stage1 = q_stage1,
  shift = b[["c"]], meanlog = b[["m"]], sdlog = b[["s"]], model_mean = b[["c"]] + E,
  stage1_mean = bmi_summary$mean, stage1_sd = bmi_summary$sd, n = bmi_summary$n)

## 4. Model profiles per unit mean effect ---------------------------------------------------
# Each profile is scaled so that the exact mean effect is 1 (dE/dc = 1, dE/dm = E,
# dE/ds^2 = E / 2). 'combined' and 'vqtl_log' take half of the mean effect from each part.
# 'vqtl_raw' is a raw-scale effect on mean and SD at constant coefficient of variation
# (Y -> (1 + e) Y), whose profile is q / mu: on this baseline (c = 0) it equals the
# multiplicative profile exactly.
additive <- W[, "additive"]; log_additive <- W[, "log_additive"] / E; variance <- W[, "sigma2"] / (E / 2)
model_mean <- b[["c"]] + E
P <- cbind(additive = additive, log_additive = log_additive, combined = 0.5 * additive + 0.5 * log_additive,
           vqtl_raw = q_model / model_mean, vqtl_log = 0.5 * log_additive + 0.5 * variance,
           variance_only = variance)
# Spans of the fitted models ("reject model b when the truth is profile a"). span{1, X}
# equals span{1, q}: the combined model and the raw-scale mean + variance model coincide.
spans <- list(additive = W[, "additive", drop = FALSE], log_additive = W[, "log_additive", drop = FALSE],
              combined = W[, c("additive", "log_additive")], vqtl_raw = cbind(1, q_model - model_mean),
              vqtl_log = W[, c("log_additive", "sigma2")], full = W)

Sigma0_inv_half <- sym_pow(Sigma0, -0.5)
to_metric <- list(whitened = function(v) as.numeric(Sigma0_inv_half %*% v), tau = function(v) as.numeric(v))
angle <- function(a, b) {
  cs <- sum(a * b) / sqrt(sum(a^2) * sum(b^2))
  c(cos = cs, deg = acos(max(-1, min(1, cs))) * 180 / pi)
}
angle_to_span <- function(a, B) {                           # angle from a to its projection on span(B)
  B <- as.matrix(B); p <- B %*% qr.solve(B, a)
  acos(min(1, sqrt(sum(p^2) / sum(a^2)))) * 180 / pi
}

profiles <- rbindlist(lapply(colnames(P), function(pn) {
  whitened_norm <- sqrt(sum((Sigma0_inv_half %*% P[, pn])^2))
  data.table(baseline = "lognormal", profile = pn, tau = u, beta_unit_mean = P[, pn],
             beta_unit_whitened_norm = P[, pn] / whitened_norm, beta_unit_tau_norm = P[, pn] / sqrt(sum(P[, pn]^2)))
}))

## 5. Angles between profiles and power to separate them ------------------------------------
# A variant with profile a and whitened length |Z| = ||sqrt(d) Sigma0^{-1/2} beta|| rejects
# model span b with non-centrality |Z|^2 sin^2(angle(a, span b)). We report the |Z| needed for
# 80% power for a targeted 1-df test (alpha 0.05 and 5e-8) and for the omnibus lack-of-fit test
# (chi-square on 45 - rank(b) df, since the whitened noise is the identity), and the sample-size
# multiple relative to a just genome-wide-significant variant (|Z| = 5.45).
ncp_for_power <- function(df, alpha) uniroot(function(l) pchisq(qchisq(1 - alpha, df), df, ncp = l,
                                                                lower.tail = FALSE) - 0.8, c(0, 1e5))$root
pair_rows <- span_rows <- separation_rows <- list()
for (metric in names(to_metric)) {
  mp <- to_metric[[metric]]
  cand <- colnames(P)
  for (i in seq_along(cand)) for (j in seq_along(cand)) if (i < j) {
    a <- angle(mp(P[, cand[i]]), mp(P[, cand[j]]))
    pair_rows[[length(pair_rows) + 1L]] <- data.table(baseline = "lognormal", metric = metric,
      profile_a = cand[i], profile_b = cand[j], cosine = a[["cos"]], angle_deg = a[["deg"]])
  }
  for (pn in cand) for (sn in names(spans)) {
    theta <- angle_to_span(mp(P[, pn]), apply(spans[[sn]], 2, mp))
    span_rows[[length(span_rows) + 1L]] <- data.table(baseline = "lognormal", metric = metric,
      truth_profile = pn, model_span = sn, span_rank = ncol(spans[[sn]]), angle_to_span_deg = theta)
    if (metric == "whitened" && theta > 1e-4) {
      s2 <- sin(theta * pi / 180)^2
      df_lof <- n_tau - ncol(spans[[sn]])
      z1 <- sqrt(ncp_for_power(1, 0.05) / s2); z1_gws <- sqrt(ncp_for_power(1, 5e-8) / s2)
      z_lof <- sqrt(ncp_for_power(df_lof, 0.05) / s2)
      separation_rows[[length(separation_rows) + 1L]] <- data.table(baseline = "lognormal", truth_profile = pn,
        model_span = sn, angle_to_span_deg = theta, z_required_targeted_1df_a05 = z1,
        z_required_targeted_1df_a5e8 = z1_gws, z_required_lof_a05 = z_lof, lof_df = df_lof,
        n_multiple_vs_z545_targeted_a05 = (z1 / 5.45)^2, n_multiple_vs_z545_lof_a05 = (z_lof / 5.45)^2)
    }
  }
}

write_result(baseline_table, OUT, "baseline")
write_result(weight_gates, OUT, "weight_gates")
write_result(profiles, OUT, "profiles")
write_result(rbindlist(pair_rows), OUT, "angles_pairwise")
write_result(rbindlist(span_rows), OUT, "angles_to_span")
write_result(rbindlist(separation_rows), OUT, "separation")
