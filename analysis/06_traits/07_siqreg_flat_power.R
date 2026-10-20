#!/usr/bin/env Rscript
# Flat Box-Cox power of BMI from non-genetic covariates and from the polygenic score
#
# What: SIQReg (Huang, Costantino and Dahl) chooses one Box-Cox scale per phenotype from 43
#   non-genetic covariates (age, age^2, sex, PC1-PC40) and reports lambda = -0.21 for BMI; our
#   flat power from genetic effects is about -1.5. Here the estimator is held fixed and only the
#   predictor changes. Each covariate, and the leave-one-chromosome-out BMI polygenic score
#   (PRS), gets an unconditional quantile-effect profile on raw BMI in the Stage 1 BMI sample
#   (07_siqreg_covariate_rif.py), and its flat power is the Box-Cox lambda at which that profile
#   has no linear tilt: the root nearest zero of
#       T(lambda) = sum_t w2_t Q(u_t)^(lambda - 1) b(u_t),   w2 = trapezoid weights x (2u - 1),
#   the same solver as 02_flat_powers.R. Per predictor we report the flat power with a
#   delta-method SE and the interval where |T| <= 1.96 SE(T); per SIQReg covariate set
#   (k = 3, 13, 23, 33, 43) the median power and two pooled powers from the precision-weighted
#   mean profile. The 639 BMI lead profiles are solved the same way for comparison.
# Paper: Discussion (lambda for BMI from covariates and from the polygenic score); Methods
#   "Change of scale"
# Input: the PRS analysis sample of the BMI polygenic-score analysis (BMI, age, sex, PC1-PC10,
#   PRS01-PRS22), UK Biobank genetic PCs 1-40 and the identifier linker; aggregate
#   results/01_data/{lead_profiles,baseline_quantiles} and results/06_traits/scale/boxcox_lambda.csv
# Output: results/06_traits/siqreg/{per_predictor,covariate_set,comparison,profiles_marginal,
#   profiles_adjusted,pooled_profiles}.csv
# Runs: on a cluster with UK Biobank data (about one hour; the RIF step dominates)

repo <- Sys.getenv("GMOQ_REPO", system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE))
source(file.path(repo, "R", "load.R"))

cfg <- read_config()
work <- file.path(cfg$bp$private_root, "06_traits", "siqreg")   # participant-level files stay here
u <- tau_grid(cfg); n_tau <- length(u)
cw <- trap_weights(u); w2 <- cw * (2 * u - 1)
lambda_grid <- round(seq(-500, 500) / 100, 2)
PCS40 <- paste0("PC", 1:40); PCS10 <- paste0("PC", 1:10); PRS <- sprintf("PRS%02d", 1:22)
COVS43 <- c("age", "age2", "sex", PCS40)
SETS <- list(k3 = COVS43[1:3], k13 = COVS43[1:13], k23 = COVS43[1:23], k33 = COVS43[1:33], k43 = COVS43)
# SIQReg's own BMI lambda by covariate set (their Supplementary Fig. 6, read off the figure) and
# with a polygenic score as the predictor (their Fig. 3b).
SIQREG_BMI <- c(k3 = 1.361, k13 = 0.157, k23 = -0.102, k33 = -0.180, k43 = -0.211)
SIQREG_BMI_PGS <- -1.3247

## 1. Analysis sample (participant level; written only to the private work directory)
# The Stage 1 BMI sample with its leave-one-chromosome-out PRS comes from the polygenic-score
# analysis; PC11-PC40 are added from the UK Biobank PC file through the identifier linker.
dir.create(file.path(work, "prep"), recursive = TRUE, showWarnings = FALSE)
A <- as.data.table(readRDS(cfg$prs$analysis_sample))     # FID IID bmi log_bmi age sex PC1-10 PRS01-22
linker <- unique(fread(cmd = paste("zcat <", shQuote(cfg$siqreg$linker)), select = c("ieu", "eid"),
                       colClasses = list(character = "ieu")), by = "eid")
pcs <- merge(fread(cfg$siqreg$pcs40, select = c("eid", PCS40)), linker, by = "eid")[, eid := NULL]
setnames(pcs, "ieu", "IID")
S <- merge(A, pcs, by = "IID", sort = FALSE, suffixes = c("", "_raw"))
S <- S[match(A$IID, S$IID)]
stopifnot(nrow(S) == nrow(A), !anyNA(S$PC40), nrow(S) == cfg$prs$n)
# The linked PC1-PC10 must equal the Stage 1 covariates exactly (checks the linkage).
stopifnot(max(vapply(PCS10, function(p) max(abs(S[[p]] - S[[paste0(p, "_raw")]])), 0)) < 1e-6)
keep <- c("bmi", "age", "sex", PCS40, PRS)
stopifnot(!anyNA(S[, keep, with = FALSE]))
fwrite(S[, keep, with = FALSE], file.path(work, "prep", "rif_input.csv.gz"), compress = "gzip")
rm(A, S, pcs, linker); invisible(gc())

## 2. Quantile-effect profiles of every predictor (Stage 1 package RIF)
helper <- file.path(repo, "analysis", "06_traits", "07_siqreg_covariate_rif.py")
status <- system2(cfg$package$stage1$python, c(helper, "--run", work, "--tau-from", cfg$grid$tau_from,
                  "--tau-to", cfg$grid$tau_to, "--n-tau", cfg$grid$n_tau),
                  stdout = file.path(work, "rif.log"), stderr = file.path(work, "rif.log"))
if (status != 0L) stop("RIF helper failed; see ", file.path(work, "rif.log"))
profiles <- fread(file.path(work, "rif", "profiles.csv"))
sigma0_all <- fread(file.path(work, "rif", "sigma0.csv"))
nuisance <- fread(file.path(work, "rif", "nuisance.csv"))
stopifnot(nrow(nuisance) == n_tau, max(abs(nuisance$tau - u)) < 1e-9)
Q_raw <- nuisance$q_raw
baseline <- as.data.table(read_result("01_data/baseline_quantiles", "baseline_quantiles"))[trait == "bmi"]
q_baseline <- as.numeric(baseline$q[match(round(u, 3), round(as.numeric(baseline$tau), 3))])
stopifnot(max(abs(Q_raw - q_baseline)) < 0.01, min(diff(Q_raw)) > 0, min(Q_raw) > 0)

## 3. The flat-power solver
# A[l, ] = w2 * Q^(lambda_l - 1), so T over the whole lambda grid is one matrix product.
tilt_matrix <- function(q, ww) outer(seq_along(lambda_grid), seq_along(q), function(i, t) ww[t] * q[t]^(lambda_grid[i] - 1))
A_full <- tilt_matrix(Q_raw, w2)
root_nearest_zero <- function(t) {
  k <- which(diff(sign(t)) != 0)
  if (!length(k)) return(if (t[1] > 0) -Inf else Inf)
  k <- k[which.min(abs(lambda_grid[k]))]
  lambda_grid[k] - t[k] * (lambda_grid[k + 1] - lambda_grid[k]) / (t[k + 1] - t[k])
}
flat_lambda <- function(b, A = A_full) root_nearest_zero(as.numeric(A %*% b))
# Delta-method SE at the root, and the run of lambda (containing the root, or nearest the
# smallest |z|) where the tilt is not significant at 5%.
flat_lambda_interval <- function(b, S0, d, lambda) {
  tilt <- as.numeric(A_full %*% b)
  tilt_z <- tilt / sqrt(rowSums((A_full %*% S0) * A_full) / d)
  se <- NA_real_
  if (is.finite(lambda)) {
    a <- w2 * Q_raw^(lambda - 1)
    dT <- sum(a * log(Q_raw) * b)                          # dT / dlambda at the root
    if (is.finite(dT) && dT != 0) se <- sqrt(as.numeric(t(a) %*% S0 %*% a) / d) / abs(dT)
  }
  ok <- is.finite(tilt_z) & abs(tilt_z) <= qnorm(0.975)
  if (!any(ok)) return(c(se, NA_real_, NA_real_))
  runs <- rle(ok); ends <- cumsum(runs$lengths); starts <- ends - runs$lengths + 1L
  i <- which(runs$values)
  anchor <- if (is.finite(lambda)) lambda else lambda_grid[which.min(abs(tilt_z))]
  j <- i[which.min(vapply(i, function(k) min(abs(lambda_grid[starts[k]:ends[k]] - anchor)), 0))]
  c(se, lambda_grid[starts[j]], lambda_grid[ends[j]])
}
# The same solver restricted to SIQReg's window tau in [0.1, 0.9] (weights rescaled to it).
in_window <- u >= 0.1 - 1e-9 & u <= 0.9 + 1e-9
u_window <- (u[in_window] - min(u[in_window])) / diff(range(u[in_window]))
A_window <- tilt_matrix(Q_raw[in_window], trap_weights(u_window) * (2 * u_window - 1))

## 4. Flat power of each predictor
sigma0_of <- function(p, m) {
  x <- sigma0_all[predictor == p & mode == m]
  as.matrix(x[order(tau), paste0("s_", sprintf("%02d", 0:(n_tau - 1))), with = FALSE])
}
profile_of <- function(p, m) { x <- profiles[predictor == p & mode == m][order(tau)]; stopifnot(nrow(x) == n_tau); x }
predictors <- unique(profiles$predictor)
per_predictor <- rbindlist(lapply(predictors, function(p) rbindlist(lapply(c("marginal", "adjusted"), function(m) {
  x <- profile_of(p, m); b <- x$b; S0 <- sigma0_of(p, m); d <- x$d[1]
  mean_z <- sum(cw * b) / sqrt(as.numeric(t(cw) %*% S0 %*% cw) / d)   # strength of the predictor
  lambda <- flat_lambda(b)
  ci <- flat_lambda_interval(b, S0, d, lambda)
  data.table(predictor = p, mode = m, kind = if (grepl("^PRS", p)) "genetic_prs" else "covariate",
             n = x$n[1], d_resid = d, sd_predictor = x$sd_predictor[1],
             mean_effect_z = mean_z, strong = abs(mean_z) > 10, lambda_flat = lambda, lambda_flat_se = ci[1], lambda_flat_lo95 = ci[2],
             lambda_flat_hi95 = ci[3], lambda_flat_tau_0.1_0.9 = flat_lambda(b[in_window], A_window),
             tilt_z_at_raw = sum(w2 * b) / sqrt(as.numeric(t(w2) %*% S0 %*% w2) / d),
             beta_mean_raw = sum(cw * b))
}))))

## 5. Per covariate set: per-predictor summaries and pooled powers
# Pooled profile: d-weighted mean of the profiles, each aligned on the sign of its mean effect.
pooled_profile <- function(ps, m) {
  B <- t(vapply(ps, function(p) profile_of(p, m)$b, numeric(n_tau)))
  d <- vapply(ps, function(p) profile_of(p, m)$d[1], 0)
  sgn <- sign(as.numeric(B %*% cw)); sgn[sgn == 0] <- 1
  as.numeric(colSums(d * sgn * B) / sum(d))
}
# Implied natural scale g(u) = int dQ / pbar and the Box-Cox power that best matches it (as in
# 02_flat_powers.R); undefined unless pbar stays positive, which it need not for weak PCs.
implied_scale_lambda <- function(pbar) {
  if (any(!is.finite(pbar)) || min(pbar) <= 0) return(c(NA_real_, NA_real_))
  g <- cumsum(c(0, diff(Q_raw) / ((pbar[-1] + pbar[-n_tau]) / 2)))
  fit_r <- vapply(lambda_grid, function(l) suppressWarnings(cor(boxcox(Q_raw, l), g)), 0)
  if (all(is.na(fit_r))) return(c(NA_real_, NA_real_))
  c(lambda_grid[which.max(fit_r)], max(fit_r, na.rm = TRUE))
}
set_row <- function(label, ps, m) {
  x <- per_predictor[predictor %in% ps & mode == m]
  clipped <- pmin(pmax(x$lambda_flat, -5), 5)            # no root: counted at the search edge
  w <- abs(x$mean_effect_z)
  pbar <- pooled_profile(ps, m)
  implied <- implied_scale_lambda(pbar)
  data.table(set = label, mode = m, n_predictors = length(ps),
             n_strong = sum(x$strong), max_abs_mean_effect_z = max(abs(x$mean_effect_z)),
             lambda_median = median(clipped), lambda_q25 = quantile(clipped, .25, names = FALSE),
             lambda_q75 = quantile(clipped, .75, names = FALSE),
             lambda_strength_weighted_mean = sum(w * clipped) / sum(w),
             lambda_median_strong = if (any(x$strong)) median(clipped[x$strong]) else NA_real_,
             lambda_pooled_tilt = flat_lambda(pbar),
             lambda_pooled_implied_scale = implied[1], implied_scale_r = implied[2],
             frac_no_root = mean(!is.finite(x$lambda_flat)))
}
covariate_set <- rbindlist(c(
  lapply(names(SETS), function(k) set_row(k, SETS[[k]], "marginal")),
  lapply(names(SETS), function(k) set_row(k, SETS[[k]], "adjusted")),
  list(set_row("prs_loco", PRS, "marginal"), set_row("prs_loco", PRS, "adjusted"))))
covariate_set[, siqreg_bmi_lambda0 := SIQREG_BMI[set]]
covariate_set[set == "prs_loco", siqreg_bmi_lambda0 := SIQREG_BMI_PGS]

## 6. The genetic comparator: the 639 BMI lead profiles, same solver
lead_profiles <- as.data.table(read_result("01_data/lead_profiles", "lead_profiles"))
leads <- lead_set_matrix(lead_profiles, "bmi_A_regenie", u)
sigma0_leads <- p05_sigma0(leads$source, u)
B_lead <- leads$B; d_lead <- leads$meta$d_resid; chr_lead <- as.character(leads$meta$chr)
sgn_lead <- sign(as.numeric(B_lead %*% cw)); sgn_lead[sgn_lead == 0] <- 1
z_lead <- sgn_lead * as.numeric(B_lead %*% cw) /
  sqrt(vapply(seq_along(d_lead), function(k) functional_cov(cw, cw, sigma0_leads[[chr_lead[k]]], d_lead[k]), 0))
lambda_lead <- apply(B_lead, 1, flat_lambda)
pbar_lead <- as.numeric(colSums(d_lead * sgn_lead * B_lead) / sum(d_lead))

summarise_powers <- function(label, lambda, mean_z) {
  clipped <- pmin(pmax(lambda, -5), 5); strong <- abs(mean_z) > 10
  data.table(predictor_class = label, n = length(lambda), n_strong = sum(strong),
             lambda_median = median(clipped), lambda_q25 = quantile(clipped, .25, names = FALSE),
             lambda_q75 = quantile(clipped, .75, names = FALSE),
             lambda_median_strong = if (any(strong)) median(clipped[strong]) else NA_real_)
}
k43_m <- per_predictor[predictor %in% COVS43 & mode == "marginal"]
k43_a <- per_predictor[predictor %in% COVS43 & mode == "adjusted"]
prs_m <- per_predictor[kind == "genetic_prs" & mode == "marginal"]
prs_a <- per_predictor[kind == "genetic_prs" & mode == "adjusted"]
comparison <- rbindlist(list(
  summarise_powers("covariates_k43_marginal", k43_m$lambda_flat, k43_m$mean_effect_z),
  summarise_powers("covariates_k43_adjusted", k43_a$lambda_flat, k43_a$mean_effect_z),
  summarise_powers("prs_loco_marginal", prs_m$lambda_flat, prs_m$mean_effect_z),
  summarise_powers("prs_loco_adjusted", prs_a$lambda_flat, prs_a$mean_effect_z),
  summarise_powers("p05_bmi_leads", lambda_lead, z_lead)))
pick <- function(s, m, col) covariate_set[set == s & mode == m][[col]]
comparison[, lambda_pooled_implied_scale := c(
  pick("k43", "marginal", "lambda_pooled_implied_scale"), pick("k43", "adjusted", "lambda_pooled_implied_scale"),
  pick("prs_loco", "marginal", "lambda_pooled_implied_scale"), pick("prs_loco", "adjusted", "lambda_pooled_implied_scale"),
  implied_scale_lambda(pbar_lead)[1])]
comparison[, lambda_pooled_tilt := c(
  pick("k43", "marginal", "lambda_pooled_tilt"), pick("k43", "adjusted", "lambda_pooled_tilt"),
  pick("prs_loco", "marginal", "lambda_pooled_tilt"), pick("prs_loco", "adjusted", "lambda_pooled_tilt"),
  flat_lambda(pbar_lead))]
print(comparison)
# Anchors: the lead-based pooled power should reproduce lambda_flat_pooled of BMI in
# boxcox_lambda.csv, and the adjusted PRS power should be close to it.
bmi_pooled <- as.data.table(read_result("06_traits/scale", "boxcox_lambda"))[trait == "bmi"]$lambda_flat_pooled
message("lead pooled power ", implied_scale_lambda(pbar_lead)[1], " vs ", bmi_pooled,
        "; adjusted PRS median ", median(prs_a$lambda_flat))

## 7. Write
out <- "06_traits/siqreg"
write_result(per_predictor, out, "per_predictor")
write_result(covariate_set, out, "covariate_set")
write_result(comparison, out, "comparison")
write_result(profiles[mode == "marginal", .(predictor, tau, b, se, d, n)], out, "profiles_marginal")
write_result(profiles[mode == "adjusted", .(predictor, tau, b, se, d, n)], out, "profiles_adjusted")
write_result(data.table(tau = u, q_raw = Q_raw, q_raw_P03 = q_baseline, fhat = nuisance$fhat,
                        pbar_k43_marginal = pooled_profile(COVS43, "marginal"),
                        pbar_prs_adjusted = pooled_profile(PRS, "adjusted"), pbar_p05_leads = pbar_lead),
             out, "pooled_profiles")
