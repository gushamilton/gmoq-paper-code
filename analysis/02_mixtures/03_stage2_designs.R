#!/usr/bin/env Rscript
# Weight-design simulation, part 2: fit the mixture and compare the weight designs
#
# What: for each simulated replicate of part 1 (03_simulate_stage1.py), fits the two-component
#   skew-normal mixture to the phenotype exactly as in the empirical analyses and estimates the
#   Stage 2 parameters by ordinary least squares with the pooled covariance under four weight
#   designs that share the same Stage 1 output:
#     par    the model design: every column -dF/dtheta divided by the fitted mixture density;
#     hyb    the hybrid design: location columns as 'par'; the membership (gamma) column
#            divided by the Stage 1 kernel density instead;
#     match  every column divided by the Stage 1 kernel density (shown for completeness);
#     oracle the true generating components, every column divided by the kernel density (what
#            the Stage 1 profile equals to first order).
#   It then summarises bias, empirical s.d. against the pooled standard error, coverage, null
#   calibration and the density of the fitted mixture against the truth in the gap.
#   gamma = logit P(lower component); the paper reports upper-component membership, -gamma.
# Paper: Supplementary Note 2 (Supplementary Note Table 1); Supplementary Fig. S7b
# Input: part 1 output per replicate; scenario definitions (03_weight_design_scenarios.json)
# Output: results/02_mixtures/weight_design/sim/{operating_characteristics, null_calibration,
#   mixing_matrix, density_by_tau, weights_by_tau, fit_summary}.csv
# Runs: on a cluster, called by 03_weight_design_simulation.sh:
#   --task --scenario S --rep-start a --rep-end b --in DIR --out DIR   (replicates a to b - 1)
#   --summary --in DIR                                                  (all replicates)

repo <- Sys.getenv("GMOQ_REPO", system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE))
source(file.path(repo, "R", "load.R"))
cfg <- read_config()
suppressPackageStartupMessages(library(gmoq))   # the paper used a development version of this package (Stage 2 commit 9fee67d)
args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag) { i <- match(flag, args); if (is.na(i)) NULL else args[[i + 1L]] }
`%||%` <- function(a, b) if (is.null(a)) b else a
scenarios <- jsonlite::read_json(file.path(repo, "analysis", "02_mixtures", "03_weight_design_scenarios.json"),
                                 simplifyVector = FALSE)
DESIGNS <- c("par", "match", "hyb", "oracle")

## 1. Distribution of a (possibly multi-part) true component ------------------------------
parts_pdf <- function(x, parts) {
  w <- vapply(parts, `[[`, 0, "w"); w <- w / sum(w)
  Reduce(`+`, lapply(seq_along(parts), function(j)
    w[j] * sn::dsn(x, xi = parts[[j]]$xi, omega = parts[[j]]$omega, alpha = parts[[j]]$alpha)))
}
parts_cdf <- function(x, parts) {
  w <- vapply(parts, `[[`, 0, "w"); w <- w / sum(w)
  Reduce(`+`, lapply(seq_along(parts), function(j)
    w[j] * sn::psn(x, xi = parts[[j]]$xi, omega = parts[[j]]$omega, alpha = parts[[j]]$alpha)))
}

## 2. Fitted baselines and their model-design weights ------------------------------------
# Mixture: the same EM fit as the empirical analyses, seeded with the replicate's seed.
fit_mixture_weights <- function(y, q_tau, taus, seed) {
  set.seed(seed)
  fit <- mixsmsn::smsn.mix(y = y, g = 2, family = "Skew.normal", nu = 1,
                           get.init = TRUE, iter.max = 500, error = 1e-4, calc.im = FALSE)
  ord <- order(as.numeric(fit$mu))
  pii <- as.numeric(fit$pii)[ord]; mu <- as.numeric(fit$mu)[ord]
  sigma2 <- as.numeric(fit$sigma2)[ord]; shape <- as.numeric(fit$shape)[ord]
  if (length(pii) != 2L || any(!is.finite(c(pii, mu, sigma2, shape))) || any(pii <= 0) || any(sigma2 <= 0))
    stop("invalid skew-normal mixture fit")
  components <- data.table(component = c("lower", "upper"), weight = pii, xi = mu,
                           omega = sqrt(sigma2), alpha = shape)
  mf <- mixture_functions(components)
  W <- gmoq::make_weights_generic(taus = taus, q_tau = q_tau, dist_cdf = mf$cdf,
                                           dist_pdf = mf$pdf, params = mf$p, grad_funcs = mf$grads)
  list(W = W, f = pmax(mf$pdf(q_tau, mf$p), 1e-12), Fq = mf$cdf(q_tau, mf$p), p = mf$p,
       iter = if (is.null(fit$iter)) NA_real_ else as.numeric(fit$iter))
}

# Shifted lognormal (BMI-like control): shift = 1st percentile - 0.001, moment estimates of
# meanlog and sdlog; columns additive (shift) and multiplicative (meanlog).
fit_lognormal_weights <- function(y, q_tau, taus) {
  shift0 <- as.numeric(stats::quantile(y, probs = 0.01, type = 8, names = FALSE)) - 1e-3
  ys <- y - shift0; ys <- ys[ys > 0]
  p <- list(shift = shift0, meanlog = mean(log(ys)), sdlog = stats::sd(log(ys)))
  pdf <- function(x, p) { o <- numeric(length(x)); ok <- x > p$shift
    o[ok] <- dlnorm(x[ok] - p$shift, p$meanlog, p$sdlog); o }
  cdf <- function(x, p) { o <- numeric(length(x)); ok <- x > p$shift
    o[ok] <- plnorm(x[ok] - p$shift, p$meanlog, p$sdlog); o }
  grad_shift <- function(x, p) -pdf(x, p)
  grad_meanlog <- function(x, p) { o <- numeric(length(x)); ok <- x > p$shift
    z <- (log(x[ok] - p$shift) - p$meanlog) / p$sdlog; o[ok] <- -dnorm(z) / p$sdlog; o }
  W <- gmoq::make_weights_generic(taus = taus, q_tau = q_tau, dist_cdf = cdf, dist_pdf = pdf,
    params = p, grad_funcs = list(additive = grad_shift, multiplicative = grad_meanlog))
  list(W = W, f = pmax(pdf(q_tau, p), 1e-12), Fq = cdf(q_tau, p), p = p)
}

## 3. Stage 2 by OLS with the pooled covariance -------------------------------------------
# theta = A b, A = (W'W)^-1 W'; Var(theta) = A Sigma0 A' / d; lack of fit Q = d r' Sigma0^-1 r.
ols_fit <- function(W, B, Sigma0, d) {
  A <- solve(crossprod(W), t(W)); V_unit <- A %*% Sigma0 %*% t(A)
  TH <- B %*% t(A); SE <- sqrt(outer(1 / d, diag(V_unit)))
  R <- B - TH %*% t(W); Si <- solve(Sigma0)
  list(A = A, Vu = V_unit, TH = TH, SE = SE, lof = d * rowSums((R %*% Si) * R))
}

read_replicate <- function(dir) {
  nu <- jsonlite::read_json(file.path(dir, "nuis.json"), simplifyVector = TRUE)
  s1 <- fread(file.path(dir, "stage1.tsv"))
  S0 <- as.matrix(fread(file.path(dir, "sigma0.tsv"), header = FALSE)); S0 <- (S0 + t(S0)) / 2
  y <- readBin(file.path(dir, "y.f64"), "double", n = nu$n, size = 8, endian = "little")
  list(nu = nu, s1 = s1, S0 = S0, y = y)
}

## 4. One replicate -----------------------------------------------------------------------
process_replicate <- function(sc, variants, dir, rep) {
  r <- read_replicate(dir)
  taus <- r$nu$taus; q <- r$nu$q_tau; f_kde <- r$nu$fhat
  beta_cols <- grep("^beta_tau_", names(r$s1), value = TRUE)
  B <- as.matrix(r$s1[, ..beta_cols]); d <- r$s1$d_resid
  ids <- r$s1$snp; vi <- match(ids, vapply(variants, `[[`, "", "id"))
  t0 <- proc.time()[["elapsed"]]
  if (sc$family == "mixture") {
    mf <- tryCatch(fit_mixture_weights(r$y, q, taus, r$nu$seed), error = function(e) e)
    if (inherits(mf, "error")) return(list(fail = data.table(rep = rep, error = conditionMessage(mf))))
    W_par <- mf$W
    pl <- sc$p_lower; lo <- sc$components$lower; up <- sc$components$upper
    FL <- parts_cdf(q, lo); FU <- parts_cdf(q, up)
    f_true <- pl * parts_pdf(q, lo) + (1 - pl) * parts_pdf(q, up)
    # Oracle: the true membership column -p(1 - p)(F_L - F_U) and component shares, over f_kde.
    W_oracle <- cbind(gamma = -pl * (1 - pl) * (FL - FU) / f_kde, beta_left = pl * parts_pdf(q, lo) / f_kde,
                      beta_right = (1 - pl) * parts_pdf(q, up) / f_kde)
    W_match <- W_par * (mf$f / f_kde)
    W_hyb <- W_par; W_hyb[, "gamma"] <- W_match[, "gamma"]
    models <- list(skew = list(par = W_par, match = W_match, hyb = W_hyb, oracle = W_oracle))
    truth_of <- function(v) c(gamma = v$gamma %||% 0, beta_left = v$beta_left %||% 0,
                              beta_right = v$beta_right %||% 0)
    fitted <- data.table(rep = rep, p1 = mf$p$p1, xi1 = mf$p$xi1, omega1 = mf$p$omega1, alpha1 = mf$p$alpha1,
                         xi2 = mf$p$xi2, omega2 = mf$p$omega2, alpha2 = mf$p$alpha2, iter = mf$iter,
                         lower_share = if (is.null(r$nu$lower_share)) NA_real_ else r$nu$lower_share,
                         max_abs_WL_plus_WR_minus_1_hyb = max(abs(W_hyb[, 2] + W_hyb[, 3] - 1)))
  } else {
    mf <- fit_lognormal_weights(r$y, q, taus)
    ln <- sc$lognormal
    f_true <- dlnorm(q - ln$shift, ln$meanlog, ln$sdlog)
    W_par <- mf$W; W_match <- W_par * (mf$f / f_kde); W_hyb <- W_par   # no membership column
    W_oracle <- cbind(additive = f_true / f_kde, multiplicative = f_true * (q - ln$shift) / f_kde)
    models <- list(additive = list(par = W_par[, 1, drop = FALSE], match = W_match[, 1, drop = FALSE],
                                   hyb = W_hyb[, 1, drop = FALSE], oracle = W_oracle[, 1, drop = FALSE]),
                   both = list(par = W_par, match = W_match, hyb = W_hyb, oracle = W_oracle))
    truth_of <- function(v) c(additive = v$additive %||% 0, multiplicative = v$multiplicative %||% 0)
    fitted <- data.table(rep = rep, shift = mf$p$shift, meanlog = mf$p$meanlog, sdlog = mf$p$sdlog,
                         max_abs_W_additive_minus_1_par = max(abs(W_par[, "additive"] - 1)),
                         max_abs_hyb_minus_par = max(abs(W_hyb - W_par)))
  }
  fit_seconds <- proc.time()[["elapsed"]] - t0

  est <- list(); mixing <- list()
  for (mn in names(models)) {
    Ws <- models[[mn]]
    fits <- lapply(Ws, ols_fit, B = B, Sigma0 = r$S0, d = d)
    pars <- colnames(Ws$par)
    for (dn in names(fits)) {
      f <- fits[[dn]]
      for (i in seq_along(ids)) {
        v <- variants[[vi[i]]]; tr <- truth_of(v)[pars]
        # A one-column additive model has a defined truth only without a multiplicative effect.
        if (mn == "additive" && (v$multiplicative %||% 0) != 0) tr[] <- NA_real_
        est[[length(est) + 1L]] <- data.table(rep = rep, model = mn, design = dn, snp = ids[i],
          role = v$role, param = pars, truth = as.numeric(tr), est = f$TH[i, ], se = f$SE[i, ],
          lof = f$lof[i], joint_T = as.numeric(d[i] * t(f$TH[i, ]) %*% solve(f$Vu, f$TH[i, ])))
      }
      # Mixing matrix M = A W_oracle: expected estimate of each parameter per unit true effect.
      mixing[[length(mixing) + 1L]] <- data.table(rep = rep, model = mn, design = dn, estimated = pars,
                                                  f$A %*% Ws$oracle, unit_se = sqrt(diag(f$Vu)))
    }
  }
  dens <- data.table(rep = rep, tau = round(taus, 3), q_tau = q, f_kde = f_kde, f_fit = mf$f, f_true = f_true,
                     F_fit_minus_tau = mf$Fq - taus)
  weights <- data.table(rep = rep, tau = round(taus, 3))
  last <- models[[length(models)]]
  for (dn in DESIGNS) for (j in colnames(last[[dn]])) weights[[paste0("W_", dn, "_", j)]] <- last[[dn]][, j]
  fitted[, `:=`(fit_seconds = fit_seconds, n = r$nu$n, bandwidth = r$nu$bandwidth,
                stage1_seconds = r$nu$stage1_seconds, direct_check = r$nu$direct_check_max_rel_err)]
  list(est = rbindlist(est), mix = rbindlist(mixing, fill = TRUE), dens = dens, w = weights, fit = fitted)
}

## 5. Task mode: a batch of replicates ----------------------------------------------------
if ("--task" %in% args) {
  scn <- getarg("--scenario"); a <- as.integer(getarg("--rep-start")); b <- as.integer(getarg("--rep-end"))
  in_dir <- getarg("--in"); out_dir <- getarg("--out"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  sc <- scenarios$scenarios[[scn]]; variants <- scenarios$snp_sets[[sc$snp_set]]
  res <- list()
  for (rep in seq.int(a, b - 1L)) {
    dir <- file.path(in_dir, scn, sprintf("rep_%04d", rep))
    if (!file.exists(file.path(dir, "READY"))) stop("missing Stage 1 replicate ", dir)
    res[[length(res) + 1L]] <- process_replicate(sc, variants, dir, rep)
    unlink(dir, recursive = TRUE)
  }
  out <- list(scenario = scn, est = rbindlist(lapply(res, `[[`, "est")),
              mix = rbindlist(lapply(res, `[[`, "mix")), dens = rbindlist(lapply(res, `[[`, "dens")),
              w = rbindlist(lapply(res, `[[`, "w")), fit = rbindlist(lapply(res, `[[`, "fit"), fill = TRUE),
              fail = rbindlist(lapply(res, `[[`, "fail")))
  for (nm in c("est", "mix", "dens", "w", "fit")) if (nrow(out[[nm]])) out[[nm]][, scenario := scn]
  saveRDS(out, file.path(out_dir, sprintf("%s_%04d.rds", scn, a)))
  quit(save = "no")
}

## 6. Summary mode: operating characteristics over all replicates -------------------------
if ("--summary" %in% args) {
  L <- lapply(list.files(getarg("--in"), pattern = "\\.rds$", full.names = TRUE), readRDS)
  E <- rbindlist(lapply(L, `[[`, "est")); MX <- rbindlist(lapply(L, `[[`, "mix"), fill = TRUE)
  D <- rbindlist(lapply(L, `[[`, "dens")); FT <- rbindlist(lapply(L, `[[`, "fit"), fill = TRUE)
  WT <- rbindlist(lapply(L, `[[`, "w"), fill = TRUE)
  FL <- rbindlist(lapply(L, function(x) if (nrow(x$fail)) x$fail[, scenario := x$scenario] else NULL), fill = TRUE)
  out_dir <- "02_mixtures/weight_design/sim"

  # Per scenario, model, design, variant and parameter with a known truth.
  E[, z := est / se]
  E[, covered := abs(est - truth) <= qnorm(0.975) * se]
  S <- E[!is.na(truth), .(n_reps = .N, truth = truth[1], mean_est = mean(est), bias = mean(est - truth),
          rel_bias = if (truth[1] != 0) mean(est - truth) / truth[1] else NA_real_,
          mcse_bias = sd(est) / sqrt(.N), emp_sd = sd(est), mean_se = mean(se),
          sd_over_se = sd(est) / mean(se), coverage = mean(covered),
          mcse_coverage = sqrt(mean(covered) * (1 - mean(covered)) / .N),
          rmse = sqrt(mean((est - truth)^2)), reject_0.05 = mean(abs(z) > qnorm(0.975)),
          mean_lof = mean(lof), median_lof = median(lof)),
         by = .(scenario, model, design, role, snp, param)]
  write_result(S[order(scenario, model, role, param, design)], out_dir, "operating_characteristics")

  # Null variants pooled: z calibration, type I error, inflation; joint Wald test (chi2_k).
  NUL <- E[role == "null", .(n = .N, mean_z = mean(z), sd_z = sd(z), type1_0.05 = mean(abs(z) > qnorm(0.975)),
            type1_0.001 = mean(abs(z) > qnorm(1 - 5e-4)), lambda = median(z^2) / qchisq(0.5, 1),
            coverage = mean(covered)), by = .(scenario, model, design, param)]
  JN <- E[role == "null", .(joint_T = joint_T[1], k = .N), by = .(scenario, model, design, rep, snp)]
  JN <- JN[, .(n = .N, joint_type1_0.05 = mean(pchisq(joint_T, k, lower.tail = FALSE) < 0.05),
               joint_lambda = median(joint_T) / qchisq(0.5, k[1])), by = .(scenario, model, design)]
  write_result(merge(NUL, JN, by = c("scenario", "model", "design"), all.x = TRUE), out_dir, "null_calibration")

  write_result(MX[, c(lapply(.SD, mean), n_reps = .N), by = .(scenario, model, design, estimated),
                  .SDcols = setdiff(names(MX), c("rep", "scenario", "model", "design", "estimated"))],
               out_dir, "mixing_matrix")
  write_result(D[, .(q_tau = mean(q_tau), f_kde = mean(f_kde), f_fit = mean(f_fit), f_true = mean(f_true),
                     ratio_fit_over_true = mean(f_fit / f_true), ratio_kde_over_true = mean(f_kde / f_true),
                     ratio_fit_over_kde = mean(f_fit / f_kde), sd_ratio_fit_over_kde = sd(f_fit / f_kde),
                     F_fit_minus_tau = mean(F_fit_minus_tau), n_reps = .N), by = .(scenario, tau)],
               out_dir, "density_by_tau")
  wcols <- grep("^W_", names(WT), value = TRUE)
  write_result(WT[, lapply(.SD, mean), by = .(scenario, tau), .SDcols = wcols], out_dir, "weights_by_tau")
  fnum <- setdiff(names(FT)[vapply(FT, is.numeric, logical(1))], "rep")
  write_result(FT[, c(list(n_reps = .N), lapply(.SD, function(v) mean(v, na.rm = TRUE))), by = scenario,
                  .SDcols = fnum], out_dir, "fit_summary")
  if (nrow(FL)) write_result(FL, out_dir, "fit_failures")
}
