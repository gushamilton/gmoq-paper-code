#!/usr/bin/env Rscript
# Do the BMI model weights recover known effects? (exact check, no simulation noise)
#
# What: validates the shifted log-normal weights used for the BMI models. For a single
#   variant with a known additive and/or multiplicative effect on a BMI-like (and a
#   weight-like) log-normal phenotype, it computes the exact population Stage 1 quantile
#   effects by numerical integration, fits the baseline to the resulting quantiles exactly as
#   the real analysis does, applies the Stage 2 weights of the joint additive + multiplicative
#   model, and reports the expected estimate against the truth. Stage 1 divides by a kernel
#   density estimate, so the check is run both with the true density and with the kernel
#   estimate's expectation. The weights are a first-order map: the bias must vanish at the
#   effect size of a typical lead (effect scale 0.1).
# Paper: Methods (Stage 2 weights for BMI)
# Input: none (synthetic; truth baselines are the shifted log-normal fits to the Stage 1
#   quantiles of BMI and weight)
# Output: results/05_bmi/weights_check/expected_bias.csv
# Runs: locally in about a minute

repo <- Sys.getenv("GMOQ_REPO", unset = NA)
if (is.na(repo)) repo <- system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE)
source(file.path(repo, "R", "load.R"))
cfg <- read_config()
u <- tau_grid(cfg)
OUT <- "05_bmi/weights_check"

# Truth baselines (y = shift + exp(meanlog + sdlog e + b g) + a g) and per-allele effects.
# Effect scale 1 is FTO-sized (b X ~ 0.28 kg/m2 per allele on BMI); 0.1 is a typical lead.
SCENARIOS <- list(
  bmi_like = list(baseline = list(shift = 13.3222, meanlog = 2.591093, sdlog = 0.3267288),
                  snps = list(A = c(a = 0.13, b = 0), M = c(a = 0, b = 0.02), AM = c(a = 0.13, b = 0.02))),
  weight_like = list(baseline = list(shift = 8.039295, meanlog = 4.224771, sdlog = 0.2228376),
                     snps = list(A = c(a = 0.35, b = 0), M = c(a = 0, b = 0.012), AM = c(a = 0.35, b = 0.012))))
ROLES <- c(A = "additive", M = "multiplicative", AM = "both")
EFFECT_SCALES <- c(1, 0.25, 0.1)

# Exact expected Stage 2 estimate for one variant with minor allele frequency `maf`.
# Given genotype g, F(y | g) = plnorm(y - c - a g; m + b g, s). With the marginal quantiles
# q_u, the population slope of the RIF regression is exactly
#   E beta(u) = - sum_g P(g) (g - gbar) F(q_u | g) / (var(g) fhat(q_u)),
# with fhat the density Stage 1 divides by: the marginal density itself, or its Gaussian-kernel
# expectation at the Silverman bandwidth for N = 100,000 (the kernel estimate's own bias).
expected_estimates <- function(ln, a, b, maf = 0.3, N = 1e5) {
  g <- 0:2; pg <- dbinom(g, 2, maf); gbar <- sum(pg * g); var_g <- sum(pg * (g - gbar)^2)
  F_g <- function(y, k) plnorm(y - ln$shift - a * g[k], ln$meanlog + b * g[k], ln$sdlog)
  f_g <- function(y, k) dlnorm(y - ln$shift - a * g[k], ln$meanlog + b * g[k], ln$sdlog)
  F_marg <- function(y) Reduce(`+`, lapply(1:3, function(k) pg[k] * F_g(y, k)))
  f_marg <- function(y) Reduce(`+`, lapply(1:3, function(k) pg[k] * f_g(y, k)))
  q <- vapply(u, function(t) uniroot(function(y) F_marg(y) - t,
         c(ln$shift + 1e-9, ln$shift + exp(ln$meanlog + 12 * ln$sdlog)), tol = 1e-13)$root, 0)
  # marginal SD by numerical integration, for the kernel bandwidth
  yy <- seq(ln$shift + 1e-6, ln$shift + exp(ln$meanlog + 9 * ln$sdlog), length.out = 200001)
  dy <- yy[2] - yy[1]; fy <- f_marg(yy); mu <- sum(yy * fy) * dy; sd_y <- sqrt(sum((yy - mu)^2 * fy) * dy)
  h <- 1.06 * sd_y * N^(-0.2); tt <- seq(-8, 8, by = 0.005); wt <- dnorm(tt) * 0.005
  f_true <- f_marg(q)
  f_kde <- vapply(q, function(x) sum(wt * f_marg(x - h * tt)), 0)
  dF <- -Reduce(`+`, lapply(1:3, function(k) pg[k] * (g[k] - gbar) * F_g(q, k))) / var_g
  # Stage 2: baseline fitted to the quantiles, weights dq/dtheta at the model quantiles
  fit <- fit_shifted_lognormal(q, u)
  W <- shifted_lognormal_weights(fit, u)[, c("additive", "multiplicative")]
  out <- list()
  for (density in c("exact_density", "kde_expectation")) {
    E_beta <- dF / (if (density == "exact_density") f_true else f_kde)
    for (sc in c("raw", "log")) {
      D <- if (sc == "raw") rep(1, length(u)) else 1 / q          # exact change of scale to log
      theta <- as.numeric(solve(crossprod(W * D), t(W * D)) %*% (E_beta * D))
      out[[length(out) + 1L]] <- data.table(stage1_density = density, design = "quantile_lsq", scale = sc,
        param = c("additive", "multiplicative"), truth = c(a, b), expected = theta)
    }
  }
  x <- rbindlist(out)
  x[, bias := expected - truth][, rel_bias := ifelse(truth != 0, bias / truth, NA_real_)]
  x[, `:=`(max_f_kde_over_f = max(f_kde / f_true), min_f_kde_over_f = min(f_kde / f_true),
           lsq_max_abs_q_err = fit[["max_abs_q_err"]])]
  x
}

expected_bias <- rbindlist(lapply(names(SCENARIOS), function(sc) {
  s <- SCENARIOS[[sc]]
  rbindlist(lapply(names(s$snps), function(id) rbindlist(lapply(EFFECT_SCALES, function(k)
    expected_estimates(s$baseline, k * s$snps[[id]][["a"]], k * s$snps[[id]][["b"]])[,
      `:=`(scenario = sc, snp = id, role = ROLES[[id]], effect_scale = k)]))))
}))
# At a typical lead's effect size the relative bias of every non-zero parameter is below 5%.
typical <- expected_bias[stage1_density == "kde_expectation" & truth != 0 & effect_scale == 0.1]
stopifnot(all(abs(typical$rel_bias) < 0.05))
write_result(expected_bias, OUT, "expected_bias")
