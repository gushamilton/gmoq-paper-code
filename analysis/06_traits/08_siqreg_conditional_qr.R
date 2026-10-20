#!/usr/bin/env Rscript
# SIQReg's estimator on our data: conditional quantile regression against our unconditional RIF
#
# What: 07_siqreg_flat_power.R changed the predictor (covariates against the polygenic score)
#   but kept our estimator. This script completes a 2 x 2 on one sample (Stage 1 BMI sample,
#   raw BMI) and SIQReg's nine quantile levels tau = 0.1, ..., 0.9:
#                               unconditional RIF (ours)   conditional QR (SIQReg's)
#     polygenic score (PRS)
#     covariates (age, sex ...)
#   (1) RIF: the profiles of 07_siqreg_covariate_rif.py recomputed at the nine levels, and the
#   same flat-power solver. Because SIQReg minimises a Wald statistic of quantile dependence
#   rather than solving for a root of the tilt, the two minimisers (of |tilt| / SE and of the
#   profile's dispersion) are reported beside the root. (2) Conditional QR: for each predictor
#   and each candidate lambda, BMI is Box-Cox transformed and the smoothed quantile regression
#   (conquer; Gaussian kernel, bandwidth ((p + log n) / n)^(2/5)) is refitted at the nine
#   levels; lambda-hat is the root nearest zero of the same linear tilt, bracketed on a unit
#   grid over [-5, 5] and bisected to 1e-3. Refitting is needed because, unlike an
#   unconditional profile, a conditional slope cannot be converted to another scale afterwards.
# Paper: Discussion ("Running both estimators on our data reproduced this"); Methods "Change of
#   scale"
# Input: the analysis sample written by 07_siqreg_flat_power.R (private work directory)
# Output: results/06_traits/siqreg_conditional_qr/{rif_tau9,conditional_qr,qr_lambda_scan,
#   comparison_2x2}.csv
# Runs: on a cluster with UK Biobank data (about three core-hours; predictors in parallel)

repo <- Sys.getenv("GMOQ_REPO", system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE))
source(file.path(repo, "R", "load.R"))

cfg <- read_config()
source_work <- file.path(cfg$bp$private_root, "06_traits", "siqreg")
work <- file.path(cfg$bp$private_root, "06_traits", "siqreg_conditional_qr")
prep <- file.path(source_work, "prep", "rif_input.csv.gz")
conquer_lib <- cfg$siqreg$conquer_lib
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
# Workers are fresh processes (PSOCK) with single-threaded BLAS: a forked worker inherits a
# running OpenBLAS thread pool and can deadlock.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
suppressPackageStartupMessages(library(conquer, lib.loc = conquer_lib))

TAU9 <- seq(0.1, 0.9, 0.1); n9 <- length(TAU9)
cw9 <- trap_weights(TAU9); w29 <- cw9 * (2 * TAU9 - 1)
lambda_grid <- round(seq(-500, 500) / 100, 2)
PCS40 <- paste0("PC", 1:40); PCS10 <- paste0("PC", 1:10); PRS <- sprintf("PRS%02d", 1:22)
STAGE1 <- c("age", "sex", PCS10)
COVS43 <- c("age", "age2", "sex", PCS40)
SETS <- list(k3 = COVS43[1:3], k13 = COVS43[1:13], k23 = COVS43[1:23], k33 = COVS43[1:33], k43 = COVS43)
SIQREG_BMI <- c(k3 = 1.361, k13 = 0.157, k23 = -0.102, k33 = -0.180, k43 = -0.211)   # their Supplementary Fig. 6
SIQREG_BMI_PGS <- -1.3247                                                           # their Fig. 3b

## 1. Unconditional RIF at the nine levels
dir.create(work, recursive = TRUE, showWarnings = FALSE)
status <- system2(cfg$package$stage1$python,
                  c(file.path(repo, "analysis", "06_traits", "07_siqreg_covariate_rif.py"), "--run", work,
                    "--prep", prep, "--tau-from", 0.1, "--tau-to", 0.9, "--n-tau", n9),
                  stdout = file.path(work, "rif.log"), stderr = file.path(work, "rif.log"))
if (status != 0L) stop("RIF helper failed; see ", file.path(work, "rif.log"))
profiles <- fread(file.path(work, "rif", "profiles.csv")); sigma0_all <- fread(file.path(work, "rif", "sigma0.csv"))
nuisance <- fread(file.path(work, "rif", "nuisance.csv"))
stopifnot(max(abs(nuisance$tau - TAU9)) < 1e-9, nuisance$n[1] == cfg$prs$n)
Q9 <- nuisance$q_raw

A9 <- outer(seq_along(lambda_grid), seq_len(n9), function(i, t) w29[t] * Q9[t]^(lambda_grid[i] - 1))
root_nearest_zero <- function(t) {
  k <- which(diff(sign(t)) != 0)
  if (!length(k)) return(if (t[1] > 0) -Inf else Inf)
  k <- k[which.min(abs(lambda_grid[k]))]
  lambda_grid[k] - t[k] * (lambda_grid[k + 1] - lambda_grid[k]) / (t[k + 1] - t[k])
}
# Coefficient of variation of a nine-level profile; NA when the predictor has no effect to
# speak of (sd / |mean| is then a ratio of two noise terms).
dispersion <- function(b) {
  m <- mean(abs(b)); d <- if (m > 0) stats::sd(b) / m else NA_real_
  if (is.finite(d)) d else NA_real_
}
# Minimisers on the RIF side: beta_T(lambda) = Q^(lambda - 1) beta exactly, so no refitting.
rif_minimisers <- function(b, S0, d) {
  tilt <- as.numeric(A9 %*% b)
  z <- abs(tilt) / sqrt(rowSums((A9 %*% S0) * A9) / d)
  BT <- outer(seq_along(lambda_grid), seq_len(n9), function(i, t) Q9[t]^(lambda_grid[i] - 1) * b[t])
  dsp <- apply(BT, 1, dispersion)
  c(if (any(is.finite(z))) lambda_grid[which.min(ifelse(is.finite(z), z, Inf))] else NA_real_,
    if (any(is.finite(dsp))) lambda_grid[which.min(ifelse(is.finite(dsp), dsp, Inf))] else NA_real_,
    min(z, na.rm = TRUE))
}
sigma0_of <- function(p, m) as.matrix(sigma0_all[predictor == p & mode == m][order(tau),
                                        paste0("s_", sprintf("%02d", 0:(n9 - 1))), with = FALSE])
profile_of <- function(p, m) { x <- profiles[predictor == p & mode == m][order(tau)]; stopifnot(nrow(x) == n9); x }
rif9 <- rbindlist(lapply(unique(profiles$predictor), function(p) rbindlist(lapply(c("marginal", "adjusted"), function(m) {
  x <- profile_of(p, m); b <- x$b; S0 <- sigma0_of(p, m); d <- x$d[1]
  mins <- rif_minimisers(b, S0, d)
  lambda <- root_nearest_zero(as.numeric(A9 %*% b))
  data.table(predictor = p, mode = m, kind = if (grepl("^PRS", p)) "genetic_prs" else "covariate",
             n = x$n[1], d_resid = d,
             mean_effect_z = sum(cw9 * b) / sqrt(as.numeric(t(cw9) %*% S0 %*% cw9) / d),
             lambda_flat_tau9 = lambda, has_root = is.finite(lambda),
             lambda_rif_min_tilt_z = mins[1], lambda_rif_min_dispersion = mins[2], min_abs_tilt_z = mins[3])
}))))

## 2. Conditional quantile regression (conquer), refitted per lambda
D <- fread(prep)
D[, age2 := age^2]
y <- as.numeric(D$bmi)
stopifnot(nrow(D) == cfg$prs$n, min(y) > 0)
bandwidth <- function(p, n) ((p + log(n)) / n)^(2 / 5)          # SIQReg's choice
qr_slopes <- function(X, yv, h) vapply(TAU9, function(tt) conquer(X, yv, tau = tt, kernel = "Gaussian", h = h)$coeff[2], 0)
# Box-Cox BMI standardised to unit SD: this multiplies the tilt by a positive constant, so the
# roots are unchanged.
tilt_at <- function(lam, X, h) {
  yT <- boxcox(y, lam); yT <- (yT - mean(yT)) / stats::sd(yT)
  b <- qr_slopes(X, yT, h)
  list(tilt = sum(w29 * b), b = b)
}
qr_one <- function(p, m) {
  cols <- if (m == "marginal") p else c(p, setdiff(STAGE1, p))   # the predictor is the first column
  X <- as.matrix(D[, cols, with = FALSE]); storage.mode(X) <- "double"
  h <- bandwidth(ncol(X), nrow(X))
  coarse <- -5:5
  fits <- lapply(coarse, tilt_at, X = X, h = h)
  tilt <- vapply(fits, `[[`, 0, "tilt")
  n_fits <- length(coarse)
  scan <- data.table(predictor = p, mode = m, lambda = coarse, tilt = tilt,
                     dispersion = vapply(fits, function(f) dispersion(f$b), 0))
  # Root of the tilt nearest zero: bracket on the unit grid, then bisect.
  k <- which(diff(sign(tilt)) != 0)
  lambda <- if (!length(k)) { if (tilt[1] > 0) -Inf else Inf } else {
    k <- k[which.min(abs(coarse[k] + 0.5))]
    lo <- coarse[k]; hi <- coarse[k + 1]; f_lo <- tilt[k]
    for (it in 1:12) {
      mid <- (lo + hi) / 2; f_mid <- tilt_at(mid, X, h)$tilt; n_fits <- n_fits + 1
      if (sign(f_mid) == sign(f_lo)) { lo <- mid; f_lo <- f_mid } else hi <- mid
      if (hi - lo < 1e-3) break
    }
    (lo + hi) / 2
  }
  # Secondary criterion, closer to SIQReg's objective: least dispersion of the profile across
  # the nine levels (coarse argmin refined on a 0.25 grid in the neighbouring unit intervals).
  lambda_disp <- NA_real_
  if (any(is.finite(scan$dispersion))) {
    j <- which.min(ifelse(is.finite(scan$dispersion), scan$dispersion, Inf))
    fine <- seq(coarse[max(1, j - 1)], coarse[min(length(coarse), j + 1)], by = 0.25)
    dg <- vapply(fine, function(l) dispersion(tilt_at(l, X, h)$b), 0)
    n_fits <- n_fits + length(fine)
    if (any(is.finite(dg))) lambda_disp <- fine[which.min(dg)]
  }
  list(scan = scan, row = data.table(predictor = p, mode = m, p_design = ncol(X), bandwidth = h,
    lambda_qr = lambda, lambda_qr_dispersion = lambda_disp,
    tilt_at_raw_scale = tilt[match(1, coarse)], n_conquer_fits = n_fits * n9))
}
jobs <- rbind(data.table(p = c(COVS43, PRS), m = "marginal"), data.table(p = c(COVS43, PRS), m = "adjusted"))
cl <- parallel::makeCluster(n_cores, type = "PSOCK")
parallel::clusterExport(cl, c("D", "y", "TAU9", "n9", "w29", "STAGE1", "conquer_lib", "jobs", "boxcox",
                              "dispersion", "bandwidth", "qr_slopes", "tilt_at", "qr_one"))
parallel::clusterEvalQ(cl, {
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
  suppressPackageStartupMessages({ library(data.table); library(conquer, lib.loc = conquer_lib) })
  NULL
})
qr_results <- parallel::clusterApplyLB(cl, seq_len(nrow(jobs)), function(i) qr_one(jobs$p[i], jobs$m[i]))
parallel::stopCluster(cl)
conditional_qr <- rbindlist(lapply(qr_results, `[[`, "row"))
qr_scan <- rbindlist(lapply(qr_results, `[[`, "scan"))
conditional_qr <- merge(conditional_qr, rif9[, .(predictor, mode, kind, mean_effect_z, lambda_flat_tau9)],
                        by = c("predictor", "mode"), all.x = TRUE, sort = FALSE)

## 3. The 2 x 2 (medians over the predictors of a set; no root counts at the search edge)
cell <- function(label, ps, m) {
  r <- rif9[predictor %in% ps & mode == m]; q <- conditional_qr[predictor %in% ps & mode == m]
  clip <- function(x) median(pmin(pmax(x, -5), 5))
  rs <- r[abs(mean_effect_z) > 10]
  data.table(predictors = label, mode = m, n_predictors = length(ps),
             n_strong = nrow(rs), n_with_root = sum(r$has_root), n_strong_with_root = sum(rs$has_root),
             rif_tau9_median = clip(r$lambda_flat_tau9),
             rif_min_tilt_z_median = clip(r$lambda_rif_min_tilt_z),
             rif_min_dispersion_median = clip(r$lambda_rif_min_dispersion),
             rif_min_dispersion_median_strong = if (nrow(rs)) clip(rs$lambda_rif_min_dispersion) else NA_real_,
             qr_median = if (nrow(q)) clip(q$lambda_qr) else NA_real_,
             qr_dispersion_median = if (nrow(q)) median(q$lambda_qr_dispersion, na.rm = TRUE) else NA_real_,
             siqreg_bmi_lambda0 = if (identical(label, "prs_loco")) SIQREG_BMI_PGS else unname(SIQREG_BMI[label]))
}
comparison_2x2 <- rbindlist(c(
  lapply(names(SETS), function(k) cell(k, SETS[[k]], "marginal")),
  list(cell("k3", SETS$k3, "adjusted"), cell("prs_loco", PRS, "marginal"), cell("prs_loco", PRS, "adjusted"))))
print(comparison_2x2)

## 4. Write
out <- "06_traits/siqreg_conditional_qr"
write_result(rif9, out, "rif_tau9")
write_result(conditional_qr, out, "conditional_qr")
write_result(qr_scan, out, "qr_lambda_scan")
write_result(comparison_2x2, out, "comparison_2x2")
