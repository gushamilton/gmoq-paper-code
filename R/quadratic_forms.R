# Lack of fit of a profile to a low-dimensional frame, with exact quadratic-form p-values.
#
# A frame X (quantile levels x k, e.g. X = [1, Q(u) - mu] for a shift plus a spread effect)
# describes a profile if its least-squares residual r = (I - H) beta, H = X (X'X)^-1 X', is
# pure noise. The statistic
#   Q = d r' Sigma0^-1 r
# is a quadratic form in a Gaussian vector, so its exact null law is sum_j lambda_j chi2_1
# with lambda the eigenvalues of Sigma0^-1 M Sigma0 M, M = I - H. The lambda do not depend
# on d, so they are computed once per Sigma0. Tail probabilities come from Davies' algorithm
# (CompQuadForm), with Imhof's method as a cross-check and the Satterthwaite scaled chi-square
# as an approximation.

# Eigenvalues of the null law of the lack-of-fit statistic for frame X.
# When baseline quantiles are tied, two influence-function columns are identical up to a
# constant and Sigma0 is exactly singular; the Moore-Penrose inverse Sigma0^+ is then used.
# Q = d r' Sigma0^+ r is still a quadratic form in r ~ N(0, M Sigma0 M / d), with exact law
# sum lambda_j chi2_1, lambda = eig(Sigma0^+ M Sigma0 M).
#   X: frame; Sigma0: pooled covariance; tol: relative eigenvalue cut-off for the rank.
# Value: list(lambda = non-zero eigenvalues, M = I - H, Si = (pseudo-)inverse, rank).
lof_spectrum <- function(X, Sigma0, tol = 1e-10) {
  n <- nrow(Sigma0)
  H <- X %*% solve(crossprod(X), t(X))
  M <- diag(n) - H
  S <- (Sigma0 + t(Sigma0)) / 2
  e <- eigen(S, symmetric = TRUE)
  rk <- sum(e$values > tol * max(e$values))
  Si <- if (rk == n) solve(S) else
    e$vectors[, 1:rk, drop = FALSE] %*% (t(e$vectors[, 1:rk, drop = FALSE]) / e$values[1:rk])
  lam <- Re(eigen(Si %*% M %*% S %*% M, only.values = TRUE)$values)
  lam <- lam[lam > 1e-10 * max(lam)]
  list(lambda = lam, M = M, Si = Si, rank = rk)
}

# Davies upper-tail p-value P(sum lambda_j chi2_1 > Q). Davies cannot resolve probabilities
# below its accuracy, so p-values under acc are set to acc and flagged rather than reported
# as 0 or negative.
#   Q: statistic(s); lambda: weights; acc, lim: Davies accuracy and iteration limit.
# Value: list(p, ifault = Davies error code, floored = TRUE where p was set to acc).
davies_p <- function(Q, lambda, acc = 1e-12, lim = 1e5) {
  if (!requireNamespace("CompQuadForm", quietly = TRUE)) stop("lof needs CompQuadForm")
  out <- vapply(Q, function(q) {
    r <- CompQuadForm::davies(q, lambda, lim = lim, acc = acc)
    c(p = r$Qq, ifault = r$ifault)
  }, numeric(2))
  p <- out["p", ]; ifault <- out["ifault", ]
  floored <- p < acc
  p[floored] <- acc
  list(p = p, ifault = ifault, floored = floored)
}

# Imhof upper-tail p-value of sum lambda_j chi2_1 (numerical integration; cross-check).
imhof_p <- function(Q, lambda, epsabs = 1e-14, epsrel = 1e-14) {
  vapply(Q, function(q) CompQuadForm::imhof(q, lambda, epsabs = epsabs, epsrel = epsrel,
                                            limit = 1e4)$Qq, numeric(1))
}

# Satterthwaite approximation: Q / c ~ chi2(nu) with c = sum(lambda^2) / sum(lambda) and
# nu = sum(lambda)^2 / sum(lambda^2) (matches the first two moments).
satterthwaite_p <- function(Q, lambda) {
  s1 <- sum(lambda); s2 <- sum(lambda^2)
  stats::pchisq(Q / (s2 / s1), s1^2 / s2, lower.tail = FALSE)
}

# Lack-of-fit test of each profile to frame X with the exact Davies reference.
#   B: leads x quantile levels (or one profile); Sigma0: pooled covariance; d: per lead;
#   X: frame; imhof_check: also return Imhof p-values.
# Value: data.frame with Q, p_davies, davies_ifault, davies_floored, p_satterthwaite,
#   df_nominal = n_tau - k, df_eff = sum(lambda)^2 / sum(lambda^2) (and p_imhof); the
#   weights lambda are attached as attribute "lambda".
lof_davies <- function(B, Sigma0, d, X, imhof_check = FALSE) {
  if (is.null(dim(B))) B <- matrix(B, nrow = 1L)
  sp <- lof_spectrum(X, Sigma0)
  R <- B %*% t(sp$M)
  Q <- d * rowSums((R %*% sp$Si) * R)
  dv <- davies_p(Q, sp$lambda)
  out <- data.frame(Q = Q, p_davies = dv$p, davies_ifault = dv$ifault,
                    davies_floored = dv$floored,
                    p_satterthwaite = satterthwaite_p(Q, sp$lambda),
                    df_nominal = nrow(Sigma0) - ncol(X),
                    df_eff = sum(sp$lambda)^2 / sum(sp$lambda^2))
  if (imhof_check) out$p_imhof <- imhof_p(Q, sp$lambda)
  attr(out, "lambda") <- sp$lambda
  out
}
