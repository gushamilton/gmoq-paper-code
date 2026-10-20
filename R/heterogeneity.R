# Between-locus heterogeneity in the shape of quantile effects.
#
# For each lead i take a location summary L_i (e.g. the change in mean) and a shape summary
# S_i (e.g. the change in spread), measured with known sampling variances vLL_i, vSS_i and
# covariance vLS_i. If every locus moved the distribution in the same way, S_i = r L_i for
# one common ratio r. The model allows the ratio to vary between loci,
#   S_i = rho_i L_i + e_i,   rho_i ~ (r, tau^2),
# and tau is the between-locus standard deviation of the ratio.
#   r   : iteratively reweighted least squares, corrected for measurement error in L,
#         r = sum w (L S - vLS) / sum w (L^2 - vLL),  w = 1 / v,  v = vSS - 2 r vLS + r^2 vLL;
#   Q   = sum (S - r L)^2 / v, referred to chi2(n - 1) (Cochran-type test of tau = 0);
#   tau = sqrt(max(0, (Q - (n - 1)) / sum((L^2 - vLL) / v)))  (method of moments).
# The alternative excess = "n" subtracts n instead of n - 1 and is reported as a sensitivity.

# Heterogeneity test and method-of-moments tau.
#   L, S: location and shape summaries per lead; vLL, vSS, vLS: their sampling (co)variances
#   (scalars are recycled; 1, 1, 0 for unit-variance uncorrelated z-scores);
#   excess: "n_minus_1" (primary) or "n"; max_iter, tol: iteration control for r.
# Value: list(r, Q, df = n - 1, p, tau, tau_over_r = tau / |r|, excess, converged,
#   resid_z = (S - r L) / sqrt(v), v, n).
het_test <- function(L, S, vLL = 1, vSS = 1, vLS = 0,
                     excess = c("n_minus_1", "n"), max_iter = 20L, tol = 1e-12) {
  excess <- match.arg(excess)
  n <- length(L)
  if (length(S) != n) stop("L and S differ in length")
  vLL <- rep_len(vLL, n); vSS <- rep_len(vSS, n); vLS <- rep_len(vLS, n)
  ok <- is.finite(L) & is.finite(S) & is.finite(vLL) & is.finite(vSS) & is.finite(vLS)
  if (!all(ok)) stop("het_test: non-finite input in ", sum(!ok), " lead(s)")
  if (n < 3L) stop("het_test needs at least 3 leads")
  r <- sum(L * S) / sum(L^2)
  converged <- FALSE
  for (it in seq_len(max_iter)) {
    v <- vSS - 2 * r * vLS + r^2 * vLL
    w <- 1 / v
    r_new <- sum(w * (L * S - vLS)) / sum(w * (L^2 - vLL))
    if (!is.finite(r_new)) break
    if (abs(r_new - r) < tol) { r <- r_new; converged <- TRUE; break }
    r <- r_new
  }
  v <- vSS - 2 * r * vLS + r^2 * vLL
  R <- S - r * L
  Q <- sum(R^2 / v)
  df <- n - 1L
  ex <- if (excess == "n") n else n - 1
  denom <- sum((L^2 - vLL) / v)
  tau <- sqrt(max(0, (Q - ex) / denom))
  list(r = r, Q = Q, df = df, p = stats::pchisq(Q, df, lower.tail = FALSE), tau = tau,
       tau_over_r = tau / abs(r), excess = excess, converged = converged,
       resid_z = R / sqrt(v), v = v, n = n)
}

# Nonparametric bootstrap over leads for tau and r.
#   B: number of resamples; seed: RNG seed; probs: percentile interval limits; other
#   arguments as het_test().
# Value: list(tau_ci, r_ci (percentile intervals), r_se, draws = B x 2 matrix (tau, r)).
het_boot <- function(L, S, vLL = 1, vSS = 1, vLS = 0, B = 2000L, seed = 20260923L,
                     excess = c("n_minus_1", "n"), probs = c(0.025, 0.975)) {
  excess <- match.arg(excess)
  n <- length(L)
  vLL <- rep_len(vLL, n); vSS <- rep_len(vSS, n); vLS <- rep_len(vLS, n)
  set.seed(seed)
  bt <- t(replicate(B, {
    i <- sample.int(n, n, replace = TRUE)
    f <- het_test(L[i], S[i], vLL[i], vSS[i], vLS[i], excess = excess)
    c(tau = f$tau, r = f$r)
  }))
  list(tau_ci = unname(stats::quantile(bt[, "tau"], probs)),
       r_ci = unname(stats::quantile(bt[, "r"], probs)),
       r_se = stats::sd(bt[, "r"]), draws = bt)
}

# Unit-variance sequential components of a profile: level, tilt given level, bend given both.
# The least-squares level, tilt and bend are theta = A beta, A = (X'X)^-1 X',
# X = [1, z, z^2] with z = qnorm(u), with covariance V0 / d, V0 = A Sigma0 A'. Gram-Schmidt
# of the seeds (anchor, e1, e2, e3) in the V0^-1 inner product gives G with G' V0^-1 G = I,
# and the components Z = sqrt(d) G^-1 theta have identity sampling covariance, so each one
# is its own z-score. c1 follows the anchor (the level by default), c2 is the tilt not
# explained by the level, and c3 the bend not explained by either.
#   B: leads x quantile levels; Sigma0: pooled covariance; d: per lead; u: quantile levels;
#   anchor: first Gram-Schmidt seed in (level, tilt, bend) coordinates.
# Value: list(Z = leads x 3 (z_c1, z_c2, z_c3), theta, G, V0, max_abs_var_minus_1 and
#   max_abs_offdiag (numerical check that the components are orthonormal)).
unit_se_components <- function(B, Sigma0, d, u, anchor = c(1, 0, 0), tol = 1e-8) {
  if (is.null(dim(B))) B <- matrix(B, nrow = 1L)
  z <- stats::qnorm(u)
  X <- cbind(level = 1, tilt = z, bend = z^2)
  A <- solve(crossprod(X), t(X))
  S <- (Sigma0 + t(Sigma0)) / 2
  V0 <- A %*% S %*% t(A)
  V0i <- chol2inv(chol((V0 + t(V0)) / 2))
  G <- NULL
  for (s in list(anchor, c(1, 0, 0), c(0, 1, 0), c(0, 0, 1))) {
    v <- as.numeric(s)
    if (!is.null(G)) for (j in seq_len(ncol(G)))
      v <- v - as.numeric(t(G[, j]) %*% V0i %*% s) * G[, j]
    nn <- as.numeric(t(v) %*% V0i %*% v)
    if (nn <= tol * as.numeric(t(s) %*% V0i %*% s)) next
    G <- cbind(G, v / sqrt(nn))
    if (ncol(G) == 3L) break
  }
  if (is.null(G) || ncol(G) != 3L) stop("unit_se_components: Gram-Schmidt failed")
  Th <- B %*% t(A)
  Z <- (Th %*% t(solve(G))) * sqrt(d)
  colnames(Z) <- c("z_c1", "z_c2", "z_c3")
  Vg <- solve(G) %*% V0 %*% t(solve(G))
  list(Z = Z, theta = Th, G = G, V0 = V0,
       max_abs_var_minus_1 = max(abs(diag(Vg) - 1)),
       max_abs_offdiag = max(abs(Vg[upper.tri(Vg)])))
}

# Detection floor of het_test on unit-variance components: the smallest tau on grid at which
# the test (Q > chi2_{1 - alpha}(n - 1)) rejects with probability at least power. Data are
# simulated as c1 = L + N(0, 1) and c2 = (r + N(0, tau^2)) L + N(0, 1), taking the observed
# c1 as the true level, which makes the floor slightly optimistic. The scan stops after three
# consecutive grid points reach the target power.
#   L: observed c1; r: fitted common ratio; grid: tau values; nsim: simulations per value.
# Value: list(floor (NA if not reached), qcrit, nsim, curve = data.frame(tau, power)).
het_detection_floor <- function(L, r, grid = seq(0, 0.5, by = 0.005), nsim = 1000L,
                                alpha = 0.05, power = 0.8, seed = 20260923L) {
  set.seed(seed)
  n <- length(L)
  qcrit <- stats::qchisq(1 - alpha, n - 1L)
  pw <- rep(NA_real_, length(grid))
  for (g in seq_along(grid)) {
    tt <- grid[g]
    pw[g] <- mean(replicate(nsim, {
      rho <- r + stats::rnorm(n, 0, tt)
      c2s <- rho * L + stats::rnorm(n)
      c1s <- L + stats::rnorm(n)
      het_test(c1s, c2s)$Q > qcrit
    }))
    if (pw[g] >= power && g >= 3L && all(pw[(g - 2L):g] >= power)) break
  }
  hit <- which(pw >= power)[1]
  list(floor = if (is.na(hit)) NA_real_ else grid[hit], qcrit = qcrit, nsim = nsim,
       curve = data.frame(tau = grid, power = pw)[!is.na(pw), ])
}
