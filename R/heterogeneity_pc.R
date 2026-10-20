# Scale-invariant heterogeneity test in the principal-component frame.
#
# The same model and estimator as het_test() (R/heterogeneity.R), but the location and
# shape summaries are each lead's scores on the first two principal directions of the
# whitened profiles,
#   y_i = sqrt(d_i) W beta_i,   W = whiten_frame(mean Sigma0)$to_white (one frame per trait),
# without scaling rows to unit length. Directions come from the noise-corrected, uncentred
# second moment of the whitened cloud,
#   M = Y'Y / n - mean_i C_i,   C_i = W Sigma0_chr(i) W'  (exact sampling covariance of y_i).
# To avoid scoring leads on directions fitted to their own noise, the leads are split at
# random into K folds and each fold is scored on directions estimated from the other folds
# (K = 2: each half scores the other). Each fold has its own common ratio r, because error in
# its directions adds a fold-specific constant to S, so
#   Q = sum_f Q_f ~ chi2(n - K),   tau = sqrt(max(0, (Q - (n - K)) / sum_f D_f)),
# with D_f = sum over fold f of (L^2 - vLL) / v. The measurement (co)variances are exact:
# vLL = v1' C_i v1, vSS = v2' C_i v2, vLS = v1' C_i v2.
#
# Why this is scale invariant: a change of scale (beta_T = D beta, Sigma_T = D Sigma0 D) maps
# y_T = O y with O orthogonal and C_T = O C O'. The directions rotate with the cloud, so
# scores, Q, p and tau are unchanged; the sign conventions below use only rotation-invariant
# quantities. The chi2(n - K) reference is slightly anti-conservative at K = 2 (fold Q's are
# correlated; het_pc_null() calibrates it) and much more so for K > 2, so K = 2 is used.

# Whitened cloud and exact per-lead sampling covariances.
#   B: leads x quantile levels; S0: named list of per-chromosome Sigma0 (as from
#   p05_sigma0(), optionally with "mean", the whitening frame); d: per lead; chr: per lead.
# Value: list(Y = whitened profiles, C = covariance per chromosome, key = chromosome per
#   lead, rank = whitening rank).
pc_cloud <- function(B, S0, d, chr) {
  if (is.null(dim(B))) B <- matrix(B, nrow = 1L)
  chr <- as.character(chr)
  if (!all(unique(chr) %in% names(S0))) stop("pc_cloud: Sigma0 missing for some chromosomes")
  S0bar <- if ("mean" %in% names(S0)) S0[["mean"]] else Reduce(`+`, S0[unique(chr)]) / length(unique(chr))
  wf <- whiten_frame(S0bar)
  W <- wf$to_white
  Y <- (B %*% t(W)) * sqrt(d)
  keys <- unique(chr)
  C <- lapply(stats::setNames(keys, keys), function(k) { M <- W %*% S0[[k]] %*% t(W); (M + t(M)) / 2 })
  list(Y = Y, C = C, key = chr, rank = wf$rank)
}

# First two principal directions of a cloud from its noise-corrected second moment
# M = Y'Y / n - mean C. Signs: PC1 points along the mean of the (aligned, positive) scores;
# PC2 is signed so that sum s1^2 s2 > 0, a rotation-invariant third moment.
# Value: list(v1, v2, values = eigenvalues of M).
pc_directions <- function(Y, key, C) {
  Cbar <- Reduce(`+`, Map(function(k, m) C[[k]] * m, names(table(key)), as.numeric(table(key)))) / length(key)
  M <- crossprod(Y) / nrow(Y) - Cbar
  e <- eigen((M + t(M)) / 2, symmetric = TRUE)
  v1 <- e$vectors[, 1]; v2 <- e$vectors[, 2]
  s1 <- Y %*% v1
  if (sum(s1) < 0) { v1 <- -v1; s1 <- -s1 }
  if (sum(s1^2 * (Y %*% v2)) < 0) v2 <- -v2
  list(v1 = v1, v2 = v2, values = e$values)
}

# Score the rows in test on directions estimated from the rows in train.
# Value: list(L = PC1 scores, S = PC2 scores, vLL, vSS, vLS (exact per lead), dirs).
pc_scores <- function(Y, key, C, train, test) {
  dr <- pc_directions(Y[train, , drop = FALSE], key[train], C)
  ks <- unique(key[test])
  q <- vapply(ks, function(k) c(vLL = sum(dr$v1 * (C[[k]] %*% dr$v1)), vSS = sum(dr$v2 * (C[[k]] %*% dr$v2)),
                                vLS = sum(dr$v1 * (C[[k]] %*% dr$v2))), numeric(3))
  j <- match(key[test], ks)
  list(L = as.numeric(Y[test, , drop = FALSE] %*% dr$v1), S = as.numeric(Y[test, , drop = FALSE] %*% dr$v2),
       vLL = q["vLL", j], vSS = q["vSS", j], vLS = q["vLS", j], dirs = dr)
}

# Random balanced assignment of n leads to K folds.
pc_folds <- function(n, seed, K = 2L) { set.seed(seed); sample(rep(seq_len(K), length.out = n)) }

# Cross-fitted heterogeneity fit: fold f is scored on directions from all other folds and
# het_test() is applied within each fold; Q and the tau denominator are summed over folds.
#   folds: fold label per lead; idx: optional list of K row-index vectors (the bootstrap
#   passes resampled indices).
# Value: list(tau, Q, df = n - K, p, n, K, r_fold, tau_fold, Q_fold, n_fold, converged,
#   per = per-fold fits and scores, idx).
het_pc_fit <- function(Y, key, C, folds, idx = NULL) {
  if (is.null(idx)) idx <- lapply(sort(unique(folds)), function(f) which(folds == f))
  K <- length(idx)
  per <- lapply(seq_len(K), function(f) {
    sc <- pc_scores(Y, key, C, train = unlist(idx[-f]), test = idx[[f]])
    h <- het_test(sc$L, sc$S, sc$vLL, sc$vSS, sc$vLS, excess = "n_minus_1")
    list(h = h, denom = sum((sc$L^2 - sc$vLL) / h$v), sc = sc)
  })
  g <- function(nm) vapply(per, function(x) x$h[[nm]], numeric(1))
  Q <- sum(g("Q")); n <- sum(g("n")); df <- n - K
  denom <- sum(vapply(per, `[[`, numeric(1), "denom"))
  tau <- sqrt(max(0, (Q - df) / denom))
  list(tau = tau, Q = Q, df = df, p = stats::pchisq(Q, df, lower.tail = FALSE), n = n, K = K,
       r_fold = g("r"), tau_fold = g("tau"), Q_fold = g("Q"), n_fold = g("n"),
       converged = all(vapply(per, function(x) x$h$converged, NA)), per = per, idx = idx)
}

# Scale-invariant heterogeneity test (main entry).
#   B: leads x quantile levels (aligned profiles); S0: per-chromosome Sigma0 list (+ "mean");
#   d: per lead; chr: chromosome per lead; seed: fold split seed; nboot: bootstrap resamples;
#   n_split: further random splits to report sensitivity to the split; K: number of folds.
# Bootstrap (nboot > 0): the primary interval holds the cross-fitted scores fixed, resamples
# leads within folds and refits r per fold (the analogue of het_boot()). Re-estimating the
# directions on each resample (tau_ci_reest, a sensitivity) attenuates tau further, because
# duplicated leads carry duplicated noise into the direction estimate.
# Value: list(tau, Q, df, p, n, r_fold, tau_fold, Q_fold, n_fold, converged, K, folds, cloud,
#   full_dirs (directions from all leads, for description and simulation only), seed, rank,
#   resid_z, L, S; plus tau_ci, boot, tau_ci_reest, boot_reest when nboot > 0 and split
#   (n_split x 2 matrix of tau and p) when n_split > 0).
het_test_pc <- function(B, S0, d, chr, seed = 20260923L, nboot = 0L, n_split = 0L, K = 2L,
                        probs = c(0.025, 0.975)) {
  cl <- pc_cloud(B, S0, d, chr)
  n <- nrow(cl$Y)
  if (n < 3L * K) stop("het_test_pc needs at least 3 leads per fold")
  folds <- pc_folds(n, seed, K)
  fit <- het_pc_fit(cl$Y, cl$key, cl$C, folds)
  idx <- fit$idx
  full <- pc_directions(cl$Y, cl$key, cl$C)
  out <- list(tau = fit$tau, Q = fit$Q, df = fit$df, p = fit$p, n = n, r_fold = fit$r_fold,
              tau_fold = fit$tau_fold, Q_fold = fit$Q_fold, n_fold = fit$n_fold, converged = fit$converged, K = K,
              folds = folds, cloud = cl, full_dirs = full, seed = seed, rank = cl$rank,
              resid_z = {
                rz <- numeric(n)
                for (f in seq_len(K)) rz[idx[[f]]] <- fit$per[[f]]$h$resid_z
                rz
              },
              L = { x <- numeric(n); for (f in seq_len(K)) x[idx[[f]]] <- fit$per[[f]]$sc$L; x },
              S = { x <- numeric(n); for (f in seq_len(K)) x[idx[[f]]] <- fit$per[[f]]$sc$S; x })
  if (nboot > 0L) {
    set.seed(seed + 1L)
    bt <- t(vapply(seq_len(nboot), function(b) {
      ib <- lapply(idx, function(i) sample.int(length(i), length(i), replace = TRUE))
      Qd <- vapply(seq_len(K), function(f) {
        sc <- fit$per[[f]]$sc; j <- ib[[f]]
        h <- het_test(sc$L[j], sc$S[j], sc$vLL[j], sc$vSS[j], sc$vLS[j], excess = "n_minus_1")
        c(h$Q, sum((sc$L[j]^2 - sc$vLL[j]) / h$v))
      }, numeric(2))
      Q <- sum(Qd[1, ]); c(tau = sqrt(max(0, (Q - (n - K)) / sum(Qd[2, ]))), Q = Q)
    }, numeric(2)))
    out$tau_ci <- unname(stats::quantile(bt[, "tau"], probs))
    out$boot <- bt
    set.seed(seed + 2L)
    br <- t(vapply(seq_len(nboot), function(b) {
      ib <- lapply(idx, function(i) i[sample.int(length(i), length(i), replace = TRUE)])
      f <- het_pc_fit(cl$Y, cl$key, cl$C, folds, idx = ib)
      c(tau = f$tau, Q = f$Q)
    }, numeric(2)))
    out$tau_ci_reest <- unname(stats::quantile(br[, "tau"], probs))
    out$boot_reest <- br
  }
  if (n_split > 0L) {
    sp <- t(vapply(seq_len(n_split), function(s) {
      f <- het_pc_fit(cl$Y, cl$key, cl$C, pc_folds(n, seed + 1000L + s, K))
      c(tau = f$tau, p = f$p)
    }, numeric(2)))
    out$split <- sp
  }
  out
}

# Simulation in the whitened frame for the null, the detection floor and recovery:
#   y_i = a_i (v1 + rho_i v2) + e_i,   rho_i ~ N(r0, tau^2),   e_i ~ N(0, C_i),
# with v1, v2 the full-data directions, a_i the observed full-data PC1 score (taken as the
# truth, which is slightly optimistic) and r0 the mean of the fitted fold ratios. The test is
# refitted on the same folds as the observed fit.
het_pc_sim_setup <- function(fit) {
  cl <- fit$cloud
  list(cl = cl, n = nrow(cl$Y), r = ncol(cl$Y), v1 = fit$full_dirs$v1, v2 = fit$full_dirs$v2,
       a = as.numeric(cl$Y %*% fit$full_dirs$v1), r0 = mean(fit$r_fold), R = lapply(cl$C, chol),
       by_key = split(seq_len(nrow(cl$Y)), cl$key), folds = fit$folds)
}

# One simulated data set at between-locus SD tt, refitted. Value: c(p, tau, Q).
het_pc_sim_once <- function(ss, tt) {
  rho <- ss$r0 + stats::rnorm(ss$n, 0, tt)
  Ys <- outer(ss$a, ss$v1) + outer(ss$a * rho, ss$v2)
  for (k in names(ss$by_key)) {
    i <- ss$by_key[[k]]
    Ys[i, ] <- Ys[i, ] + matrix(stats::rnorm(length(i) * ss$r), length(i)) %*% ss$R[[k]]
  }
  f <- het_pc_fit(Ys, ss$cl$key, ss$cl$C, ss$folds)
  c(p = f$p, tau = f$tau, Q = f$Q)
}

# Simulated null distribution (tau = 0) of Q for this trait's leads and folds. Each fold's Q
# is chi2(n_f - 1) given its directions, but the fold Q's are positively correlated (each
# fold's directions carry the other folds' noise), so Q is over-dispersed relative to
# chi2(n - K). The calibrated p-value matches the null mean m and SD s with a scaled
# chi-square, Q ~ c chi2(nu), c = s^2 / (2 m), nu = 2 m^2 / s^2, which keeps a chi-square tail.
#   fit: output of het_test_pc(); nsim: null simulations; alpha: level for q_crit.
# Value: list(Q_null, mean, sd, dispersion = sd / sqrt(2 df), q_crit = null 1 - alpha
#   quantile, size_chi2 = rejection rate of the chi2(n - K) reference, p_empirical, scale_c,
#   nu, p_calibrated, log10_p_calibrated).
het_pc_null <- function(fit, nsim = 1000L, alpha = 0.05, seed = 20260923L) {
  ss <- het_pc_sim_setup(fit)
  set.seed(seed)
  q0 <- replicate(nsim, het_pc_sim_once(ss, 0)[["Q"]])
  m <- mean(q0); s <- stats::sd(q0)
  list(Q_null = q0, mean = m, sd = s, dispersion = s / sqrt(2 * fit$df),
       q_crit = unname(stats::quantile(q0, 1 - alpha)), size_chi2 = mean(stats::pchisq(q0, fit$df, lower.tail = FALSE) < alpha),
       p_empirical = (1 + sum(q0 >= fit$Q)) / (nsim + 1),
       scale_c = s^2 / (2 * m), nu = 2 * m^2 / s^2,
       p_calibrated = stats::pchisq(fit$Q / (s^2 / (2 * m)), 2 * m^2 / s^2, lower.tail = FALSE),
       log10_p_calibrated = stats::pchisq(fit$Q / (s^2 / (2 * m)), 2 * m^2 / s^2, lower.tail = FALSE, log.p = TRUE) / log(10))
}

# Detection floor of het_test_pc(): the smallest tau on grid with power at least `power`,
# simulating the full statistic (direction estimation and cross-fitting included). The scan
# stops after three consecutive grid points reach the target power.
#   q_crit: critical value of Q (the calibrated one from het_pc_null()); NULL uses chi2(n - K).
# Value: list(floor, r0, nsim, q_crit, curve = data.frame(tau, power)).
het_pc_detection_floor <- function(fit, grid = seq(0, 0.5, by = 0.005), nsim = 1000L,
                                   alpha = 0.05, power = 0.8, seed = 20260923L, q_crit = NULL) {
  ss <- het_pc_sim_setup(fit)
  if (is.null(q_crit)) q_crit <- stats::qchisq(1 - alpha, fit$df)
  set.seed(seed)
  pw <- rep(NA_real_, length(grid))
  for (g in seq_along(grid)) {
    pw[g] <- mean(replicate(nsim, het_pc_sim_once(ss, grid[g])[["Q"]] > q_crit))
    if (pw[g] >= power && g >= 3L && all(pw[(g - 2L):g] >= power)) break
  }
  hit <- which(pw >= power)[1]
  list(floor = if (is.na(hit)) NA_real_ else grid[hit], r0 = ss$r0, nsim = nsim, q_crit = q_crit,
       curve = data.frame(tau = grid, power = pw)[!is.na(pw), ])
}

# Recovery of a known tau by het_test_pc() on this trait's leads. Noise in the cross-fitted
# directions attenuates tau (an estimated PC2 mixes in noise directions), more so when the
# second eigenvalue is small relative to the noise.
# Value: list(tau_true, mean, sd of the estimates, ratio = mean / tau_true, nsim).
het_pc_recovery <- function(fit, tau, nsim = 200L, seed = 20260923L) {
  ss <- het_pc_sim_setup(fit)
  set.seed(seed)
  est <- replicate(nsim, het_pc_sim_once(ss, tau)[["tau"]])
  list(tau_true = tau, mean = mean(est), sd = stats::sd(est), ratio = mean(est) / tau, nsim = nsim)
}
