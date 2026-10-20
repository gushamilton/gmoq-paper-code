# Shifted-lognormal baseline for the BMI and weight models.
#
# The models of a variant's effect on BMI (additive shift, multiplicative change, change in
# spread) are written for a shifted-lognormal baseline Y = c + exp(m + s Z), Z ~ N(0, 1),
# whose quantiles are q(u) = c + exp(m + s z_u), z_u = qnorm(u). The Stage 2 weight columns
# are the derivatives of these quantiles with respect to each model parameter.

# Fit c, m and s to the Stage 1 baseline quantiles by least squares:
#   minimise sum_t (c + exp(m + s z_t) - Q_t)^2  subject to s > 0 and c < min(Q).
# Starts: c0 = min(Q) - k (median(Q) - min(Q)) for k in 0.05, ..., 4, with m and s from the
# least-squares line of log(Q - c0) on z; each start is refined by Nelder-Mead then BFGS and
# the best fit kept. (Matching moments of log(Q - c) for a fixed shift fits BMI-like
# baselines poorly, hence the direct fit.)
#   q: baseline quantiles; u: quantile levels.
# Value: c(c, m, s, sse = residual sum of squares, conv = optim convergence code,
#   max_abs_q_err = largest absolute quantile error).
fit_shifted_lognormal <- function(q, u) {
  z <- stats::qnorm(u)
  obj <- function(p) if (p[3] <= 0 || p[1] >= min(q)) 1e10 else sum((p[1] + exp(p[2] + p[3] * z) - q)^2)
  span <- stats::median(q) - min(q)
  best <- NULL
  for (k in c(0.05, 0.25, 0.5, 1, 2, 4)) {
    c0 <- min(q) - k * span
    ab <- stats::coef(stats::lm(log(q - c0) ~ z))
    f <- stats::optim(unname(c(c0, ab[1], ab[2])), obj, control = list(maxit = 50000, reltol = 1e-15))
    f2 <- stats::optim(f$par, obj, method = "BFGS", control = list(maxit = 10000, reltol = 1e-15))
    if (f2$value < f$value) f <- f2
    if (is.null(best) || f$value < best$value) best <- f
  }
  p <- unname(best$par)
  qm <- p[1] + exp(p[2] + p[3] * z)
  c(c = p[1], m = p[2], s = p[3], sse = best$value, conv = best$convergence,
    max_abs_q_err = max(abs(qm - q)))
}

# Stage 2 weight columns dq_u / dtheta at the fitted model quantiles q_u = c + X_u,
# X_u = exp(m + s z_u): additive (d/dc) = 1, multiplicative (d/dm) = X_u, sdlog (d/ds) = X_u z_u.
#   par: output of fit_shifted_lognormal(); u: quantile levels.
# Value: n_tau x 3 matrix with columns additive, multiplicative, sdlog.
shifted_lognormal_weights <- function(par, u) {
  z <- stats::qnorm(u); X <- exp(par[["m"]] + par[["s"]] * z)
  cbind(additive = rep(1, length(u)), multiplicative = X, sdlog = X * z)
}
