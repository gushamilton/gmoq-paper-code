# Change of measurement scale applied to Stage 1 output.
#
# Quantiles are equivariant under any strictly increasing transform T: Q_T(u) = T(Q(u)).
# A quantile effect therefore converts by the chain rule,
#   beta_T(u) = T'(Q(u)) beta(u),   Sigma_T = D Sigma0 D,   D = diag(T'(Q(u))),
# with d (the per-variant scaling of the sampling covariance) unchanged. This is exact apart
# from the kernel density estimate used inside the Stage 1 influence function, so a quantile
# GWAS on one scale gives the profiles on every other scale without refitting.

# Derivative T'(Q(u)) of a change of scale at the baseline quantiles.
#   q: baseline quantiles Q(u) on the current scale;
#   transform: "identity" (1), "log" (1 / q), "boxcox" (q^(lambda - 1), the derivative of
#   (y^lambda - 1) / lambda), or "normal_score" (fhat(Q(u)) / phi(qnorm(u)), the map to a
#   normal-score reference scale on which the baseline is standard normal);
#   lambda: Box-Cox parameter; fhat: baseline density at q; u: quantile levels.
# Value: vector D, one entry per quantile level.
scale_derivative <- function(q, transform = c("identity", "log", "boxcox", "normal_score"),
                             lambda = NULL, fhat = NULL, u = NULL) {
  transform <- match.arg(transform)
  switch(transform,
    identity = rep(1, length(q)),
    log = { if (any(q <= 0)) stop("log scale needs positive quantiles"); 1 / q },
    boxcox = {
      if (is.null(lambda)) stop("boxcox needs lambda")
      if (any(q <= 0)) stop("Box-Cox needs positive quantiles")
      q^(lambda - 1)
    },
    normal_score = {
      if (is.null(fhat) || is.null(u)) stop("normal_score needs fhat and u")
      fhat / stats::dnorm(stats::qnorm(u))
    })
}

# Apply a change of scale: beta_T = beta D (column-wise) and Sigma_T = D Sigma0 D.
#   B: leads x quantile levels (or one profile); Sigma0: pooled covariance; D: derivative
#   vector from scale_derivative().
# Value: list(B = converted profiles, Sigma0 = converted covariance).
change_scale <- function(B, Sigma0, D) {
  if (is.null(dim(B))) B <- matrix(B, nrow = 1L)
  if (length(D) != ncol(B)) stop("D must have one entry per tau")
  list(B = sweep(B, 2, D, `*`), Sigma0 = (D * Sigma0) * rep(D, each = length(D)))
}

# Box-Cox transform: (y^lambda - 1) / lambda, and log(y) at lambda = 0.
boxcox <- function(y, lambda) if (abs(lambda) < 1e-12) log(y) else (y^lambda - 1) / lambda
