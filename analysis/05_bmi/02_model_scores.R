#!/usr/bin/env Rscript
# Which model best describes each BMI lead variant?
#
# What: fits seven models of a genetic effect (additive, multiplicative about zero,
#   multiplicative about a fitted origin, additive + multiplicative, log-scale location and
#   variance, all three log-normal parameters, and a model-free level/tilt/bend) to the Stage 1
#   quantile profile of each of the 639 BMI lead variants, by ordinary least squares, and scores
#   each fit by the lack-of-fit statistic Q with the full pooled covariance, plus AIC-like
#   (Q + 2k) and BIC-like (Q + k log N) penalties. Winners and their margins are descriptive.
#   The models' weights are built on a shifted log-normal baseline fitted to the Stage 1
#   quantiles; the same construction gives the weights of the height-adjusted weight trait,
#   which the genome-wide mixture analyses use.
# Paper: Supplementary Fig. S13; Supplementary Table S6; Fig. 4 (legend)
# Input: Stage 1 quantiles (results/01_data/baseline_quantiles); lead profiles and pooled
#   covariance of raw and log BMI (results/01_data/lead_profiles)
# Output: results/05_bmi/model_scores/{baseline_fit,weight_matrices,weight_gates,
#   model_definitions,lead_scores,lead_winners,win_counts,gap_summary}.csv
# Runs: locally from results/ in about a minute. The paper used a development version of
#   the G-MoQ R package (Stage 2 commit 9fee67d); gmoq 0.7.0 has the same estimators.

repo <- Sys.getenv("GMOQ_REPO", unset = NA)
if (is.na(repo)) repo <- system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE)
source(file.path(repo, "R", "load.R"))
source(file.path(repo, "analysis", "05_bmi", "bmi_common.R"))
suppressPackageStartupMessages(library(gmoq))   # G-MoQ Stage 2
cfg <- read_config()
u <- tau_grid(cfg); n_tau <- length(u); z <- qnorm(u)
N_BMI <- cfg$cohorts$bmi$n
OUT <- "05_bmi/model_scores"

baseline_q <- read_result("01_data/baseline_quantiles", "baseline_quantiles")
q_of <- function(trait) stage1_quantiles(baseline_q, trait, u)$q
profiles <- read_lead_profiles()

## 1. Fit the shifted log-normal baseline ------------------------------------------------
# q(u) = c + exp(m + s z_u), fitted to the 45 Stage 1 quantiles by least squares on the
# quantiles themselves, so that the weights are evaluated where the data are.
fits <- list(bmi = fit_shifted_lognormal(q_of("bmi"), u), weight = fit_shifted_lognormal(q_of("weight"), u))
baseline_fit <- rbindlist(lapply(names(fits), function(trait) {
  fp <- fits[[trait]]; q <- q_of(trait)
  qm <- fp[["c"]] + exp(fp[["m"]] + fp[["s"]] * z); mid <- which.min(abs(u - 0.5))
  data.table(trait = trait, baseline = "quantile_lsq", shift = fp[["c"]], meanlog = fp[["m"]], sdlog = fp[["s"]],
    q_err_tau002 = qm[1] - q[1], q_err_tau050 = qm[mid] - q[mid], q_err_tau098 = qm[n_tau] - q[n_tau],
    max_abs_q_err = max(abs(qm - q)), rms_q_err = sqrt(mean((qm - q)^2)),
    max_rel_q_err = max(abs(qm - q) / q), optim_convergence = fp[["conv"]])
}))

## 2. Check the weights against the package ---------------------------------------------
# The Stage 2 package builds weights generically as -dF/dtheta / f from a model's CDF and
# density. At the model quantiles this must reproduce the analytic dq/dtheta, with the
# additive column exactly 1.
weight_gates <- rbindlist(lapply(names(fits), function(trait) {
  fp <- fits[[trait]]; qm <- fp[["c"]] + exp(fp[["m"]] + fp[["s"]] * z)
  par <- list(shift = fp[["c"]], meanlog = fp[["m"]], sdlog = fp[["s"]])
  cdf <- function(x, p) { o <- numeric(length(x)); k <- x > p$shift; o[k] <- plnorm(x[k] - p$shift, p$meanlog, p$sdlog); o }
  pdf <- function(x, p) { o <- numeric(length(x)); k <- x > p$shift; o[k] <- dlnorm(x[k] - p$shift, p$meanlog, p$sdlog); o }
  d_shift <- function(x, p) -pdf(x, p)
  d_meanlog <- function(x, p) { zz <- (log(x - p$shift) - p$meanlog) / p$sdlog; -dnorm(zz) / p$sdlog }
  d_sdlog <- function(x, p) { zz <- (log(x - p$shift) - p$meanlog) / p$sdlog; -dnorm(zz) * zz / p$sdlog }
  W_package <- gmoq::make_weights_generic(u, qm, cdf, pdf, par,
                 list(additive = d_shift, meanlog = d_meanlog, sdlog = d_sdlog))
  W_analytic <- bmi_model_weights(fp, q_of(trait), u)$shifted_lognormal_full
  rel <- apply(abs(W_package - W_analytic), 2, max) / apply(abs(W_analytic), 2, max)
  add_dev <- max(abs(W_package[, "additive"] - 1))
  data.table(trait = trait, column = colnames(W_analytic), max_rel_diff_package_vs_analytic = rel,
             max_abs_Wadd_minus_1 = add_dev, pass = rel < 1e-8 & add_dev < 1e-10)
}))
stopifnot(all(weight_gates$pass))

## 3. Model weights on the raw and log scales -----------------------------------------------
# Log BMI uses the raw-scale models carried over by the exact change of scale
# W_log(u) = W_raw(u) / Q_raw(u) (d log y / dy at the quantiles), except the model-free
# level/tilt/bend basis, which is defined directly on each scale.
W_raw <- bmi_model_weights(fits$bmi, q_of("bmi"), u)
MODELS <- names(W_raw)
dlog_dy <- 1 / exp(q_of("log_bmi"))
W_log <- lapply(W_raw, function(W) W * dlog_dy)
W_log$three_components <- W_raw$three_components
W_by_scale <- list(raw = W_raw, log = W_log)

model_definitions <- rbindlist(lapply(names(W_by_scale), function(sc) { L <- W_by_scale[[sc]]
  data.table(scale = sc, model = names(L), columns = vapply(L, function(W) paste(colnames(W), collapse = ","), ""),
             k = vapply(L, ncol, 0L), rank = vapply(L, function(W) check_rank(W)$rank, 0L), set = "primary")
}))
weight_matrices <- rbindlist(c(
  lapply(names(W_by_scale), function(sc) rbindlist(lapply(MODELS, function(m)
    data.table(trait = "bmi", scale = sc, model = m, tau = u, as.data.table(W_by_scale[[sc]][[m]]))), fill = TRUE)),
  list(rbindlist(lapply(c("additive", "shifted", "both"), function(m)
    data.table(trait = "weight", scale = "raw", model = m, tau = u,
               as.data.table(bmi_model_weights(fits$weight, q_of("weight"), u)[[m]]))), fill = TRUE))), fill = TRUE)

## 4. Score every lead under every model ----------------------------------------------------
# OLS fit theta = (W'W)^-1 W' beta, residual r = (I - H) beta, and the lack-of-fit statistic
# Q = r' Sigma^-1 r with Sigma = Sigma0_chr / d (pooled covariance of the lead's chromosome).
# AIC-like = Q + 2 rank(W), BIC-like = Q + rank(W) log N. These are plug-in, descriptive scores
# (OLS rather than maximum likelihood). The Stage 2 package computes them; a direct
# computation checks it.
score_rows <- list()
for (sc in names(W_by_scale)) {
  leads <- read_lead_set(profiles, if (sc == "raw") "bmi_A_regenie" else "log_bmi_A_regenie", u)
  Sigma0 <- read_lead_sigma0(leads$source, u)
  meta <- leads$meta; B <- leads$B; Wl <- W_by_scale[[sc]]
  for (k in unique(as.character(meta$chr))) {
    i <- which(as.character(meta$chr) == k)
    Sk <- Sigma0[[k]]; d <- meta$d_resid[i]
    # compact pooled covariance: the upper triangle of Sigma0 row by row, scaled by 1 / d per lead
    iu <- which(upper.tri(Sk, diag = TRUE)); ord <- order(row(Sk)[iu], col(Sk)[iu])
    stage1 <- list(Q_slope = t(B[i, , drop = FALSE]), SE_tau = t(sqrt(outer(1 / d, diag(Sk)))))
    colnames(stage1$Q_slope) <- colnames(stage1$SE_tau) <- meta$snp_id[i]
    fit <- gmoq::param_gwas_multi(stage1, Wl, N = N_BMI, se_mode = "tau_cov",
             cov_data = list(cov_vec = as.numeric(Sk[iu[ord]]), offsets = rep(0, length(i)), cov_scales = 1 / d),
             residual_loss = "full_cov")
    Sigma0_inv <- solve((Sk + t(Sk)) / 2)
    for (m in names(Wl)) {
      X <- Wl[[m]]; Bi <- B[i, , drop = FALSE]
      R <- Bi - Bi %*% t(X %*% solve(crossprod(X), t(X)))
      Q_direct <- d * rowSums((R %*% Sigma0_inv) * R)
      stopifnot(max(abs(fit$Q_full[, m] - Q_direct) / pmax(Q_direct, 1e-12)) < 1e-6,
                all(fit$reason[, m] == "ok"))
      score_rows[[length(score_rows) + 1L]] <- data.table(scale = sc, snp_id = meta$snp_id[i], chr = meta$chr[i],
        bp = meta$bp[i], delta_mu_z = meta$delta_mu_z[i], model = m, rank = unname(fit$rank[m]),
        Q_full = fit$Q_full[, m], AIC = fit$AIC[, m], BIC = fit$BIC[, m], reason = fit$reason[, m], set = "primary")
    }
  }
}
lead_scores <- rbindlist(score_rows)

## 5. Winners and margins ------------------------------------------------------------------
winners <- function(scores, criterion) {
  x <- scores[, .(model, v = get(criterion)), by = .(scale, snp_id, chr, delta_mu_z)]
  setorder(x, scale, snp_id, v)
  x[, .(winner = model[1], runner_up = model[2], gap = v[2] - v[1], n_models = .N, winner_value = v[1]),
    by = .(scale, snp_id, chr, delta_mu_z)]
}
lead_winners <- rbindlist(lapply(c("AIC", "BIC"), function(cr) winners(lead_scores, cr)[, `:=`(criterion = cr, set = "primary")]))
lead_winners[, strength_bin := cut(abs(delta_mu_z), c(0, 7, 10, 20, 40, Inf), include.lowest = TRUE)]
win_counts <- lead_winners[, .(n = .N), by = .(scale, criterion, set, model = winner)]
win_counts[, pct := 100 * n / sum(n), by = .(scale, criterion, set)]
setorder(win_counts, scale, criterion, set, -n)
gap_summary <- lead_winners[, .(n_leads = .N, median_gap = median(gap), q25_gap = quantile(gap, 0.25),
  q75_gap = quantile(gap, 0.75), frac_gap_lt_2 = mean(gap < 2), frac_gap_lt_6 = mean(gap < 6),
  frac_gap_ge_10 = mean(gap >= 10)), by = .(scale, criterion, set)]

write_result(baseline_fit, OUT, "baseline_fit")
write_result(weight_matrices, OUT, "weight_matrices")
write_result(weight_gates, OUT, "weight_gates")
write_result(model_definitions, OUT, "model_definitions")
write_result(lead_scores, OUT, "lead_scores")
write_result(lead_winners, OUT, "lead_winners")
write_result(win_counts, OUT, "win_counts")
write_result(gap_summary, OUT, "gap_summary")
