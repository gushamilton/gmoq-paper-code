# Whitening and principal components of lead quantile-effect profiles.
#
# A lead's Stage 1 profile beta_i (45 quantile effects) has sampling covariance
# Sigma0 / d_i, where Sigma0 is the pooled covariance of the Stage 1 influence functions
# and d_i the residual genotype variance times the sample size. Whitening,
#   y_i = sqrt(d_i) Sigma0^(-1/2) beta_i,
# gives every lead identity sampling covariance, so noise is the same in every direction
# and principal components describe shape rather than the noise structure. Each row is then
# scaled to unit length (so loci contribute shape, not strength) and the second-moment
# matrix is decomposed without centring. The number of components is chosen by parallel
# analysis against isotropic noise. Loadings map back to the quantile scale as
# Sigma0^(1/2) v. Angles between directions are reported in the whitened metric.

# Symmetric matrix power S^p via the eigen-decomposition S = E diag(l) E'.
#   S: symmetric positive-definite matrix; p: power (e.g. -0.5 for the inverse square root).
# Value: E diag(l^p) E'.
sym_pow <- function(S, p) {
  S <- (S + t(S)) / 2
  e <- eigen(S, symmetric = TRUE)
  if (min(e$values) <= 0) stop("sym_pow: matrix is not positive definite")
  e$vectors %*% diag(e$values^p, length(e$values)) %*% t(e$vectors)
}

# Symmetric whitening of profiles: row i becomes sqrt(d_i) beta_i' Sigma0^(-1/2).
#   B: leads x quantile levels (or one profile); Sigma0: pooled covariance; d: per lead.
# Value: matrix of whitened profiles, one row per lead.
whiten <- function(B, Sigma0, d) {
  if (is.null(dim(B))) B <- matrix(B, nrow = 1L)
  Sih <- sym_pow(Sigma0, -0.5)
  (B %*% Sih) * sqrt(d)
}

# Whitening restricted to the column space of Sigma0. With eigenvectors E_r and eigenvalues
# L_r above tol * max eigenvalue, y = sqrt(d) L_r^(-1/2) E_r' beta has identity sampling
# covariance in rank(Sigma0) dimensions, and a whitened direction v maps back to the quantile
# scale as E_r L_r^(1/2) v. For full-rank Sigma0 this is the symmetric whitening rotated by
# E', so spectra, parallel analysis, quantile-scale loadings and angles are unchanged. When
# baseline quantiles are tied (a coarsely recorded trait) Sigma0 is exactly singular and the
# null directions, in which beta carries no information, are dropped.
#   Sigma0: pooled covariance; tol: relative eigenvalue cut-off.
# Value: list(to_white = r x n_tau matrix (y = to_white %*% beta), to_tau = n_tau x r matrix
#   (loading back to the quantile scale), E = retained eigenvectors (E %*% v gives a whitened
#   vector in symmetric-whitening coordinates), rank = r).
whiten_frame <- function(Sigma0, tol = 1e-10) {
  S <- (Sigma0 + t(Sigma0)) / 2
  e <- eigen(S, symmetric = TRUE)
  rk <- sum(e$values > tol * max(e$values))
  E <- e$vectors[, 1:rk, drop = FALSE]; L <- e$values[1:rk]
  list(to_white = t(E) / sqrt(L),
       to_tau = E * rep(sqrt(L), each = nrow(E)),
       E = E, rank = rk)
}

# Uncentred PCA of unit-length rows. Each row y_i is scaled to u_i = y_i / |y_i| and the
# (optionally weighted) second moment M = sum_i w_i u_i u_i' / sum_i w_i is decomposed.
#   Y: leads x dimensions; centre: subtract column means of the unit rows first;
#   w: optional weight per row.
# Value: list(values = eigenvalues of M, vectors, scores = U %*% vectors, norms = |y_i|,
#   pct = percentage of the total second moment per component).
pca_unit <- function(Y, centre = FALSE, w = NULL) {
  norms <- sqrt(rowSums(Y^2))
  U <- Y / norms
  if (centre) U <- sweep(U, 2, colMeans(U))
  if (is.null(w)) w <- rep(1, nrow(U))
  M <- crossprod(U * sqrt(w)) / sum(w)
  e <- eigen((M + t(M)) / 2, symmetric = TRUE)
  list(values = e$values, vectors = e$vectors, scores = U %*% e$vectors, norms = norms,
       pct = 100 * e$values / sum(e$values))
}

# Sequential parallel analysis for the number of shape components. For k = 0, ..., k_max the
# signal is taken to be the rank-k SVD reconstruction of the whitened profiles Y; unit
# isotropic noise is added, the same unit-length normalisation and uncentred PCA applied, and
# the observed share of component k + 1 is compared with its simulated null. The number of
# components is the first k whose (k + 1)th share does not exceed the null at level alpha.
#   Y: whitened profiles (not unit-normalised); k_max: largest k tried; draws: null draws per
#   k; seed: RNG seed; alpha: level; noise: optional function(n, p) returning a noise matrix.
# Value: list(table = one row per k with observed and null shares and the Monte Carlo p,
#   n_components, draws, seed).
parallel_analysis <- function(Y, k_max = 4L, draws = 300L, seed = 20260923L, alpha = 0.05,
                              noise = NULL) {
  set.seed(seed)
  obs <- pca_unit(Y)
  share <- obs$values / sum(obs$values)
  sv <- svd(Y)
  rows <- vector("list", k_max + 1L)
  for (k in 0:k_max) {
    Yk <- if (k == 0L) matrix(0, nrow(Y), ncol(Y)) else
      sv$u[, 1:k, drop = FALSE] %*% diag(sv$d[1:k], k, k) %*% t(sv$v[, 1:k, drop = FALSE])
    nullk <- replicate(draws, {
      E <- if (is.null(noise)) matrix(stats::rnorm(length(Y)), nrow(Y)) else noise(nrow(Y), ncol(Y))
      ev <- pca_unit(Yk + E)$values
      ev[k + 1L] / sum(ev)
    })
    rows[[k + 1L]] <- data.frame(
      k_assumed = k, component_tested = k + 1L,
      observed_pct = 100 * share[k + 1L],
      null_median_pct = 100 * stats::median(nullk),
      null_p95_pct = 100 * unname(stats::quantile(nullk, 0.95)),
      p_exceeds_null = (1 + sum(nullk >= share[k + 1L])) / (draws + 1))
  }
  tab <- do.call(rbind, rows)
  above <- tab$p_exceeds_null <= alpha
  n_comp <- if (all(above)) k_max + 1L else which(!above)[1] - 1L
  list(table = tab, n_components = n_comp, draws = draws, seed = seed)
}

# Acute angle in degrees between two directions: acos(|a'b| / (|a| |b|)). Eigenvector signs
# are arbitrary, so the angle is folded into [0, 90].
angle_deg <- function(a, b) {
  cs <- sum(a * b) / (sqrt(sum(a^2)) * sqrt(sum(b^2)))
  acos(min(1, abs(cs))) * 180 / pi
}

# Angle in the whitened metric between two quantile-scale directions a and b: the angle
# between Sigma0^(-1/2) a and Sigma0^(-1/2) b.
angle_whitened_deg <- function(a, b, Sigma0) {
  Sih <- sym_pow(Sigma0, -0.5)
  angle_deg(Sih %*% a, Sih %*% b)
}

# How well a loading v is described by a basis X (columns, e.g. [1, z] for shift + spread):
# R^2 = 1 - RSS / TSS of the least-squares fit of v on X. Loadings are directions, so the
# uncentred version (TSS = sum v^2) is primary; the centred version is returned alongside.
#   v: loading on the quantile scale; X: basis matrix (include a constant column if wanted).
# Value: c(r2_uncentred, r2_centred).
basis_fit_r2 <- function(v, X) {
  fit <- stats::lm.fit(X, v)
  rss <- sum(fit$residuals^2)
  c(r2_uncentred = 1 - rss / sum(v^2), r2_centred = 1 - rss / sum((v - mean(v))^2))
}
