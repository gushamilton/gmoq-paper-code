# Joint Wald tests of Stage 2 parameters and a rank check of the weight matrix.
#
# Stage 2 parameter estimates theta (k of them) have covariance V estimated from the pooled
# Stage 1 covariance, which itself has nu degrees of freedom (the number of samples behind
# the pooled matrix). The Wald statistic T = theta' V^-1 theta is therefore referred to
# Hotelling's T^2 distribution rather than to chi2(k):
#   T (nu - k + 1) / (k nu) ~ F(k, nu - k + 1).
# For the sample sizes in this paper the two references agree closely.

# Hotelling p-value of a Wald statistic.
#   T: statistic; k: number of parameters; nu: covariance degrees of freedom (passed as is;
#   the formula subtracts k itself).
# Value: upper-tail p-value.
hotelling_p <- function(T, k, nu) {
  if (any(nu <= k, na.rm = TRUE)) stop("hotelling_p: nu must exceed k")
  stats::pf(T * (nu - k + 1) / (k * nu), k, nu - k + 1, lower.tail = FALSE)
}

# Wald test T = theta' V^-1 theta with a Hotelling reference. A near-singular V (reciprocal
# condition number below rcond_min) or nu <= k is reported as refused, with T and p missing,
# rather than repaired.
#   theta: estimates; V: their covariance; nu: covariance degrees of freedom.
# Value: list(T, k, nu, p, p_chisq, refused = NA or the reason).
wald_test <- function(theta, V, nu, rcond_min = 1e-12) {
  theta <- as.numeric(theta); V <- as.matrix(V)
  k <- length(theta)
  if (!all(dim(V) == k)) stop("wald_test: V must be k x k")
  rc <- rcond(V)
  if (!is.finite(rc) || rc < rcond_min)
    return(list(T = NA_real_, k = k, nu = nu, p = NA_real_, refused = sprintf("rcond %.3g", rc)))
  if (nu <= k) return(list(T = NA_real_, k = k, nu = nu, p = NA_real_, refused = "nu <= k"))
  T <- as.numeric(t(theta) %*% solve(V, theta))
  list(T = T, k = k, nu = nu, p = hotelling_p(T, k, nu), p_chisq = stats::pchisq(T, k, lower.tail = FALSE),
       refused = NA_character_)
}

# Column rank of a Stage 2 weight matrix W by QR. Stage 2 is identified only if W has full
# column rank.
# Value: list(rank, full = TRUE when rank equals the number of columns).
check_rank <- function(W, tol = 1e-10) {
  qrW <- qr(W, tol = tol)
  list(rank = qrW$rank, full = qrW$rank == ncol(W))
}
