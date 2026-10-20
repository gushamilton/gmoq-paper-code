# Helpers for the polygenic-score interaction analysis and for projecting profiles on a basis.
# Heterogeneity tests themselves use het_test() (R/heterogeneity.R).

# Ordinary least squares with model-based and heteroskedasticity-robust (HC0) covariance.
#   b = (X'X)^-1 X'y,  e = y - X b,
#   V_model = (X'X)^-1 sum(e^2) / (n - p),  V_hc0 = (X'X)^-1 X' diag(e^2) X (X'X)^-1.
#   X: design matrix including the intercept column; y: outcome; XtXi: optional (X'X)^-1.
# Value: list(b, V_model, V_hc0, n).
ols_hc0 <- function(X, y, XtXi = NULL) {
  if (is.null(XtXi)) XtXi <- solve(crossprod(X))
  b <- drop(XtXi %*% crossprod(X, y))
  e <- drop(y - X %*% b)
  n <- nrow(X); p <- ncol(X)
  V_model <- XtXi * (sum(e^2) / (n - p))
  V_hc0 <- XtXi %*% crossprod(X * e) %*% XtXi
  dimnames(V_hc0) <- dimnames(V_model) <- list(colnames(X), colnames(X))
  names(b) <- colnames(X)
  list(b = b, V_model = V_model, V_hc0 = V_hc0, n = n)
}

# Rank-based inverse normal transform: qnorm((rank(y) - 0.5) / n), ties averaged.
rank_inverse_normal <- function(y) stats::qnorm((rank(y, ties.method = "average") - 0.5) / length(y))

# Least-squares projection of profiles on a basis H with exact sampling covariance:
#   theta_i = A beta_i,  A = (H'H)^-1 H',  Var(theta_i) = A Sigma0_i A' / d_i.
#   B: leads x quantile levels; H: quantile levels x k basis; Sigma0: one matrix, or a list
#   with one matrix per lead; d: per lead.
# Value: list(theta = leads x k, V = leads x k x k array of covariances, A).
project_profiles <- function(B, H, Sigma0, d) {
  if (is.null(dim(B))) B <- matrix(B, nrow = 1L)
  A <- solve(crossprod(H), t(H))
  Th <- B %*% t(A)
  k <- ncol(H)
  V <- array(NA_real_, c(nrow(B), k, k))
  for (i in seq_len(nrow(B))) {
    S0 <- if (is.list(Sigma0)) Sigma0[[i]] else Sigma0
    V[i, , ] <- (A %*% S0 %*% t(A)) / d[i]
  }
  list(theta = Th, V = V, A = A)
}
