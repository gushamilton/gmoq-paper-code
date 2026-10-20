# Linear functionals of a quantile-effect profile and their exact variances.
#
# Stage 1 gives, for each variant, effects beta(u_t) on the population quantiles at 45 levels
# u_1 = 0.02, ..., u_45 = 0.98, with sampling covariance Sigma0 / d (Sigma0 pooled per trait
# and chromosome, d per variant). Any summary that is linear in the profile, f = w' beta, has
#   Var(w' beta) = w' Sigma0 w / d,   Cov(a' beta, b' beta) = a' Sigma0 b / d,
# exactly, so no resampling over quantile levels is needed. The summaries used in the paper
# are trapezoid integrals over u of beta(u) times a weight function:
#   delta_mu      = int beta(u) du                      (change in the mean; equals the
#                                                       change in the first L-moment)
#   delta_lambda2 = int beta(u) (2u - 1) du             (second L-moment: spread)
#   delta_lambda3 = int beta(u) (6u^2 - 6u + 1) du      (third L-moment: skew)
#   delta_var     = 2 int beta(u) (Q(u) - mu) du        (change in the variance)
# The grid covers [0.02, 0.98], so these integrals cover 96% of the distribution; delta_mu is
# reported as computed, not rescaled by 1 / 0.96.

# Trapezoid-rule weights for integrating over the quantile levels u:
#   int f(u) du ~ sum_t c_t f(u_t),  c_1 = (u_2 - u_1) / 2,  c_n = (u_n - u_{n-1}) / 2,
#   c_t = (u_{t+1} - u_{t-1}) / 2 otherwise.
# On the 45-level grid sum(c) = 0.96.
trap_weights <- function(u) {
  n <- length(u)
  if (n < 2L) stop("trap_weights needs at least two grid points")
  d <- diff(u)
  w <- numeric(n)
  w[1] <- d[1] / 2
  w[n] <- d[n - 1] / 2
  if (n > 2L) w[2:(n - 1)] <- (d[-1] + d[-length(d)]) / 2
  w
}

# Weight vectors w (one entry per quantile level) for the named functionals f = w' beta.
#   u: quantile levels; q: baseline quantiles Q(u) on the analysed scale (optional; needed
#   for delta_var and skew_beyond_frame).
# Value: named list of weight vectors:
#   delta_mu, delta_lambda2, delta_lambda3: the L-moment integrals above (c = trapezoid weights);
#   mean_grid: the plain average over the grid, 1 / 45 each (the mean of the 45 effects);
#   delta_var: 2 c (Q - mu_b), with mu_b = sum(c Q) / sum(c) the baseline mean on the grid;
#   skew_beyond_frame: the third-L-moment weights with the shift + spread frame
#     X = [1, Q - mu_b] projected out, (I - H) w_lambda3, H = X (X'X)^-1 X'. This is skew
#     that a change of location and scale cannot produce;
#   cf_level, cf_tilt, cf_bend: rows of (X3'X3)^-1 X3' with X3 = [1, z, z^2 - 1],
#     z = qnorm(u): least-squares level, tilt and bend of the profile in normal scores.
functional_weights <- function(u, q = NULL) {
  cw <- trap_weights(u)
  z <- stats::qnorm(u)
  W <- list(
    delta_mu = cw,
    delta_lambda2 = cw * (2 * u - 1),
    delta_lambda3 = cw * (6 * u^2 - 6 * u + 1),
    mean_grid = rep(1 / length(u), length(u))
  )
  if (!is.null(q)) {
    if (length(q) != length(u)) stop("q and u differ in length")
    mu_b <- sum(cw * q) / sum(cw)
    W$delta_var <- 2 * cw * (q - mu_b)
    X <- cbind(1, q - mu_b)
    H <- X %*% solve(crossprod(X), t(X))
    W$skew_beyond_frame <- as.numeric((diag(length(u)) - H) %*% W$delta_lambda3)
  }
  X3 <- cbind(1, z, z^2 - 1)
  A3 <- solve(crossprod(X3), t(X3))
  W$cf_level <- A3[1, ]; W$cf_tilt <- A3[2, ]; W$cf_bend <- A3[3, ]
  W
}

# Estimates, standard errors and covariance of a set of functionals for each variant.
#   B: variants x quantile levels (or one profile); Sigma0: pooled covariance; d: one value
#   per variant (or one shared value); W: named list of weight vectors or an n_tau x m matrix.
# Value: list(estimate = B W (variants x m), se = sqrt(diag(W' Sigma0 W) / d_i),
#   cov_unit = W' Sigma0 W (covariance for d = 1), d, cov = function(i) giving the m x m
#   covariance for variant i).
functional_estimate <- function(B, Sigma0, d, W) {
  if (is.null(dim(B))) B <- matrix(B, nrow = 1L)
  Wm <- if (is.list(W)) do.call(cbind, W) else W
  if (is.null(colnames(Wm))) colnames(Wm) <- paste0("f", seq_len(ncol(Wm)))
  if (nrow(Wm) != ncol(B)) stop("weights and profiles differ in tau length")
  if (length(d) == 1L) d <- rep(d, nrow(B))
  if (length(d) != nrow(B)) stop("d must have one value per variant")
  est <- B %*% Wm
  V0 <- t(Wm) %*% Sigma0 %*% Wm
  se <- sqrt(outer(1 / d, diag(V0)))
  colnames(se) <- colnames(Wm)
  list(estimate = est, se = se, cov_unit = V0, d = d,
       cov = function(i) V0 / d[i])
}

# Exact covariance of two functionals a' beta and b' beta: a' Sigma0 b / d (vector over d).
functional_cov <- function(a, b, Sigma0, d) as.numeric(t(a) %*% Sigma0 %*% b) / d
