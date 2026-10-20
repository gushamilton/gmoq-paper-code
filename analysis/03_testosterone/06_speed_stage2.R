#!/usr/bin/env Rscript
# Cost of G-MoQ and of direct mixture maximum likelihood, part 2: Stage 2
#
# What: times Stage 2 alone, given the Stage 1 output of each simulated panel written by part 1
#   (06_speed_profile.py), in one warm R session on one core: building the mixture weight
#   matrix, param_gwas() with the benchmark's arguments (OLS, pooled covariance, saddlepoint
#   lack-of-fit p-values), the same without the per-variant lack-of-fit calibration, and the
#   bare projection theta = A b with its pooled standard errors. Each is timed for 100 and
#   100,000 variants, and for 1,200,000 variants at N = 400,000 (larger panels repeat the 100
#   simulated profiles; the cost does not depend on their values). Simulated data only.
# Paper: Fig. 2c; Supplementary Fig. S5d
# Input: <run_dir>/stage2_input_n<N>_rep1.json and the benchmark's Stage 2 output for checking
# Output: <run_dir>/stage2_timings_<label>.csv
# Runs: Rscript 06_speed_stage2.R <run_dir> <R library holding gmoq> <label>
#   (called by 06_speed.sh), in about 10 minutes

a <- commandArgs(trailingOnly = TRUE)
run <- a[1]; lib <- a[2]; label <- a[3]
.libPaths(unique(c(lib, .libPaths())))
elapsed <- function(expr) { t <- proc.time()[["elapsed"]]; force(expr); proc.time()[["elapsed"]] - t }
t_lib <- elapsed(suppressPackageStartupMessages(library(gmoq, lib.loc = lib)))
version <- as.character(packageVersion("gmoq"))
rows <- list()
rec <- function(n, component, M, seconds, reps = 1L, note = "")
  rows[[length(rows) + 1L]] <<- data.frame(label = label, pkg_version = version, n = n, component = component,
                                           M = M, seconds = seconds, reps = reps, note = note)
rec(0, "r_library_load_warm_process", NA, t_lib)

# Repeat the 100 simulated profiles to make a panel of M variants.
tile <- function(x, M) {
  idx <- rep(seq_len(ncol(x$Q_slope)), ceiling(M / ncol(x$Q_slope)))[seq_len(M)]
  list(stage = list(taus = x$taus, q_tau = x$q_tau, Q_slope = x$Q_slope[, idx, drop = FALSE],
                    SE_tau = x$SE_tau[, idx, drop = FALSE]),
       n_eff = x$n_eff[idx], cov_scales = x$cov_scales[idx])
}
# Stage 2 with the benchmark's arguments.
fit_benchmark <- function(s, W, cov_vec, n, q_method = "saddlepoint") {
  M <- ncol(s$stage$Q_slope)
  gmoq::param_gwas(stage1 = s$stage, transform = "custom_W", transform_args = list(W = W),
    se_mode = "tau_cov", cov_data = list(cov_vec = cov_vec, offsets = rep(0, M), cov_scales = s$cov_scales),
    tau_cov_estimator = "ols", n_threads = 1L, return_cov = "none", q_method = q_method, nu = n,
    nu_source = "panel_pooled", n_blocks = NULL, n_eff = s$n_eff, unknown_nu = "error")
}
# The projection itself: theta = A b with A = (W'W)^-1 W', Var = A Sigma0 A' / d.
bare_projection <- function(s, W, Sigma) {
  A <- solve(crossprod(W), t(W)); th <- A %*% s$stage$Q_slope
  v <- diag(A %*% Sigma %*% t(A)); list(th, sqrt(outer(v, s$cov_scales)))
}

files <- sort(list.files(run, pattern = "^stage2_input_n[0-9]+_rep1\\.json$", full.names = TRUE))
n_max <- max(as.numeric(sub(".*_n([0-9]+)_rep1\\.json$", "\\1", files)))
for (f in files) {
  t <- elapsed(x <- jsonlite::fromJSON(f, simplifyVector = TRUE)); n <- as.numeric(x$n)
  rec(n, "r_read_json_M100", 100, t)
  x$Q_slope <- as.matrix(x$beta); x$SE_tau <- as.matrix(x$se_tau)
  b <- x$baseline; taus <- as.numeric(x$taus); q <- as.numeric(x$q_tau)
  make_W <- function() gmoq::make_weights_normal_mixture(taus, q, as.numeric(b$p1), as.numeric(b$mu1),
    as.numeric(b$sd1), as.numeric(b$mu2), as.numeric(b$sd2), include_membership = TRUE)
  W <- make_W(); R <- 200L; rec(n, "r_make_W", NA, elapsed(for (i in seq_len(R)) make_W()) / R, R)
  cov_vec <- as.numeric(x$cov_vec); Sigma <- matrix(0, 45, 45); Sigma[upper.tri(Sigma, diag = TRUE)] <- cov_vec
  Sigma <- Sigma + t(Sigma) - diag(diag(Sigma))
  x$n_eff <- as.numeric(x$n_eff); x$cov_scales <- as.numeric(x$cov_scales)
  for (M in c(100L, 100000L, if (n == n_max) 1200000L)) {
    s <- tile(x, M); R <- if (M == 100L) 20L else 1L
    if (M == 100L) {   # the benchmark's own Stage 2 estimates must be reproduced
      ref <- fit_benchmark(s, W, cov_vec, n)
      bo <- jsonlite::fromJSON(sub("stage2_input", "stage2_output", f), simplifyVector = FALSE)
      est <- do.call(cbind, lapply(bo$rows, function(r) as.numeric(unlist(r$estimates))))
      rec(n, "check_max_abs_diff_vs_child", M, max(abs(est - ref$params)), 1L, "not a time: estimate agreement")
    }
    rec(n, "r_param_gwas_bench_args", M, elapsed(for (i in seq_len(R)) fit_benchmark(s, W, cov_vec, n)) / R, R,
        "pooled tau_cov OLS, saddlepoint Q p-values")
    rec(n, "r_param_gwas_chisq_naive", M, elapsed(for (i in seq_len(R)) fit_benchmark(s, W, cov_vec, n, "chisq_naive")) / R, R,
        "no per-SNP Q calibration")
    Rb <- if (M == 100L) 200L else 3L
    rec(n, "r_bare_projection_plus_se", M, elapsed(for (i in seq_len(Rb)) bare_projection(s, W, Sigma)) / Rb, Rb,
        "A %*% B and pooled SEs")
    rm(s); gc(FALSE)
  }
}
utils::write.csv(do.call(rbind, rows), file.path(run, paste0("stage2_timings_", label, ".csv")), row.names = FALSE)
