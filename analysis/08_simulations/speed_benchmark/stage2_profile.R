#!/usr/bin/env Rscript
# aux_speed_profile, part 2 (R): cost of G-MoQ Stage 2 proper, given Stage 1 output.
# Reads the Stage 2 input JSONs written by speed_profile.py (the same payload the frozen
# benchmark's --panel-stage2 child reads), and times, in one warm R session:
#   W construction (make_weights_normal_mixture), param_gwas() with the benchmark's exact
#   arguments (pooled tau_cov, OLS, saddlepoint Q p-values), the same with q_method =
#   "chisq_naive" (no per-SNP Q calibration), and the bare projection A %*% B plus pooled
#   delta-method SEs, for M = 100, 100,000 and (N = 400k input only) 1,200,000 SNPs.
# Larger M tiles the 100 benchmark columns (Stage 2 cost does not depend on the values).
# Usage: Rscript stage2_profile.R <run_dir> <r_lib> <label> [<r_package>]
#   r_package: fungwasStage2 (the paper pins, default) or gmoq (released package).
a <- commandArgs(trailingOnly = TRUE)
run <- a[1]; lib <- a[2]; label <- a[3]
pkg <- if (length(a) >= 4) a[4] else "fungwasStage2"
.libPaths(unique(c(lib, .libPaths())))
t0 <- proc.time()[["elapsed"]]
suppressPackageStartupMessages(library(pkg, lib.loc = lib, character.only = TRUE))
t_lib <- proc.time()[["elapsed"]] - t0
ver <- as.character(packageVersion(pkg))
pkg_fun <- function(name) getExportedValue(pkg, name)
el <- function(expr) { t <- proc.time()[["elapsed"]]; force(expr); proc.time()[["elapsed"]] - t }
rows <- list()
rec <- function(n, component, M, seconds, reps = 1L, note = "")
  rows[[length(rows) + 1L]] <<- data.frame(label = label, pkg_version = ver, n = n, component = component,
                                           M = M, seconds = seconds, reps = reps, note = note)
rec(0, "r_library_load_warm_process", NA, t_lib)

tile <- function(x, M) {
  k <- ceiling(M / ncol(x$Q_slope)); idx <- rep(seq_len(ncol(x$Q_slope)), k)[seq_len(M)]
  list(stage = list(taus = x$taus, q_tau = x$q_tau, Q_slope = x$Q_slope[, idx, drop = FALSE],
                    SE_tau = x$SE_tau[, idx, drop = FALSE]),
       n_eff = x$n_eff[idx], cov_scales = x$cov_scales[idx])
}
fit_bench <- function(s, W, cv, n, qm = "saddlepoint") {
  M <- ncol(s$stage$Q_slope)
  pkg_fun("param_gwas")(stage1 = s$stage, transform = "custom_W", transform_args = list(W = W),
    se_mode = "tau_cov", cov_data = list(cov_vec = cv, offsets = rep(0, M), cov_scales = s$cov_scales),
    tau_cov_estimator = "ols", n_threads = 1L, return_cov = "none", q_method = qm, nu = n,
    nu_source = "panel_pooled", n_blocks = NULL, n_eff = s$n_eff, unknown_nu = "error")
}
bare <- function(s, W, Sigma) {             # the projection itself: theta = A b, Var = A Sigma A' / d
  A <- solve(crossprod(W), t(W)); th <- A %*% s$stage$Q_slope
  v <- diag(A %*% Sigma %*% t(A)); se <- sqrt(outer(v, s$cov_scales)); list(th, se)
}

files <- sort(list.files(run, pattern = "^stage2_input_n[0-9]+_rep1\\.json$", full.names = TRUE))
for (f in files) {
  t <- el(x <- jsonlite::fromJSON(f, simplifyVector = TRUE)); n <- as.numeric(x$n)
  rec(n, "r_read_json_M100", 100, t)
  x$Q_slope <- as.matrix(x$beta); x$SE_tau <- as.matrix(x$se_tau)
  b <- x$baseline; taus <- as.numeric(x$taus); q <- as.numeric(x$q_tau)
  mk <- pkg_fun("make_weights_normal_mixture")
  # The paper's weights use the model density (config density_mode); packages that
  # offer a density option (fungwasStage2 >= 0.5.0, gmoq) are asked for it explicitly,
  # since gmoq defaults to "hybrid" when W has a membership column.
  dens <- if ("density" %in% names(formals(mk))) list(density = "model") else list()
  mkW <- function() do.call(mk, c(list(taus, q, as.numeric(b$p1), as.numeric(b$mu1),
           as.numeric(b$sd1), as.numeric(b$mu2), as.numeric(b$sd2), include_membership = TRUE), dens))
  W <- mkW(); R <- 200L; rec(n, "r_make_W", NA, el(for (i in seq_len(R)) mkW()) / R, R)
  cv <- as.numeric(x$cov_vec); Sigma <- matrix(0, 45, 45); Sigma[upper.tri(Sigma, diag = TRUE)] <- cv
  Sigma <- Sigma + t(Sigma) - diag(diag(Sigma))
  x$n_eff <- as.numeric(x$n_eff); x$cov_scales <- as.numeric(x$cov_scales)
  Ms <- c(100L, 100000L); if (n == max(as.numeric(sub(".*_n([0-9]+)_rep1\\.json$", "\\1", files)))) Ms <- c(Ms, 1200000L)
  for (M in Ms) {
    s <- tile(x, M); R <- if (M == 100L) 20L else 1L
    child_out <- sub("stage2_input", "stage2_output", f)
    if (M == 100L && file.exists(child_out)) {
      ref <- fit_bench(s, W, cv, n)   # the benchmark child's exact output must be reproduced
      bo <- jsonlite::fromJSON(child_out, simplifyVector = FALSE)
      est <- do.call(cbind, lapply(bo$rows, function(r) as.numeric(unlist(r$estimates))))
      rec(n, "check_max_abs_diff_vs_child", M, max(abs(est - ref$params)), 1L, "not a time: estimate agreement")
    }
    tq <- el(for (i in seq_len(R)) fit_bench(s, W, cv, n)) / R
    rec(n, "r_param_gwas_bench_args", M, tq, R, "pooled tau_cov OLS, saddlepoint Q p-values")
    tn <- el(for (i in seq_len(R)) fit_bench(s, W, cv, n, "chisq_naive")) / R
    rec(n, "r_param_gwas_chisq_naive", M, tn, R, "no per-SNP Q calibration")
    Rb <- if (M == 100L) 200L else 3L
    rec(n, "r_bare_projection_plus_se", M, el(for (i in seq_len(Rb)) bare(s, W, Sigma)) / Rb, Rb, "A %*% B and pooled SEs")
    rm(s); gc(FALSE)
  }
}
out <- do.call(rbind, rows)
utils::write.csv(out, file.path(run, paste0("stage2_timings_", label, ".csv")), row.names = FALSE)
print(out)
