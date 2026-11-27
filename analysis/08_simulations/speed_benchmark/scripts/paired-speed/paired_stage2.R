#!/usr/bin/env Rscript
# Public pooled-covariance Stage 2 with fitted phenotype-only mixture weights.
a <- commandArgs(trailingOnly=TRUE)
if (length(a) != 5L) stop("usage: paired_stage2.R PIN.json input.json output.json primary|reuse repetitions")
elapsed <- function() as.numeric(Sys.time())
t_entry <- elapsed()
if (!requireNamespace("jsonlite", quietly=TRUE)) stop("jsonlite is required")
pin <- jsonlite::fromJSON(a[1], simplifyVector=TRUE)
if (is.null(pin$r_lib)) stop("PIN.json must include r_lib")
.libPaths(unique(c(normalizePath(pin$r_lib, mustWork=TRUE), .libPaths())))
pkg <- if (is.null(pin$r_package)) "gmoq" else pin$r_package
t <- elapsed()
suppressPackageStartupMessages(library(pkg, lib.loc=pin$r_lib, character.only=TRUE))
library_seconds <- elapsed() - t
pkg_fun <- function(name) getExportedValue(pkg, name)
version <- as.character(packageVersion(pkg))
if (!is.null(pin$package_version) && version != pin$package_version) stop("R package version differs from pin")
t <- elapsed()
x <- jsonlite::fromJSON(a[2], simplifyVector=TRUE)
read_seconds <- elapsed() - t
stage <- list(taus=as.numeric(x$taus), q_tau=as.numeric(x$q_tau),
              Q_slope=as.matrix(x$beta), SE_tau=as.matrix(x$se_tau))
b <- x$baseline
mk <- pkg_fun("make_weights_normal_mixture")
args <- list(stage$taus, stage$q_tau, as.numeric(b$p1), as.numeric(b$mu1),
             as.numeric(b$sd1), as.numeric(b$mu2), as.numeric(b$sd2), include_membership=TRUE)
# The frozen benchmark's analytic weights use model density. State it explicitly
# for versions whose current default is hybrid density.
if ("density" %in% names(formals(mk))) args$density <- "model"
t <- elapsed()
W <- do.call(mk, args)
weight_seconds <- elapsed() - t
M <- ncol(stage$Q_slope)
source <- if (is.null(x$source)) "pooled" else x$source
pooled <- identical(source, "pooled")
fit_once <- function() pkg_fun("param_gwas")(
  stage1=stage, transform="custom_W", transform_args=list(W=W), se_mode="tau_cov",
  cov_data=list(cov_vec=as.numeric(x$cov_vec), offsets=as.numeric(x$offsets),
                cov_scales=if (pooled) as.numeric(x$cov_scales) else NULL),
  tau_cov_estimator="ols", n_threads=1L, return_cov="none", q_method="saddlepoint",
  nu=if (pooled) as.numeric(x$n) else NULL, nu_source=if (pooled) "panel_pooled" else NULL,
  n_blocks=if (pooled) NULL else 1000L, n_eff=as.numeric(x$n_eff), unknown_nu="error")
t <- elapsed()
fit <- fit_once()
fit_seconds <- elapsed() - t
if (any(!is.finite(fit$params)) || any(!is.finite(fit$SE_params)) || any(fit$SE_params<=0)) {
  stop("non-finite Stage 2 estimate or invalid SE")
}
timings <- list(library_load_seconds=library_seconds, read_json_seconds=read_seconds,
                make_weights_seconds=weight_seconds, first_param_gwas_seconds=fit_seconds)
if (a[4] == "reuse") {
  R <- as.integer(a[5])
  if (is.na(R) || R < 1L) stop("positive reuse repetition count required")
  times <- numeric(R)
  max_est_diff <- 0
  max_se_diff <- 0
  for (i in seq_len(R)) {
    t <- elapsed()
    reused <- fit_once()
    times[i] <- elapsed() - t
    max_est_diff <- max(max_est_diff, max(abs(reused$params - fit$params)))
    max_se_diff <- max(max_se_diff, max(abs(reused$SE_params - fit$SE_params)))
  }
  if (!is.finite(max_est_diff) || !is.finite(max_se_diff) || max_est_diff > 1e-12 || max_se_diff > 1e-12) {
    stop("Stage 2 reuse changed estimates or SEs")
  }
  timings$reuse_repetitions <- R
  timings$reuse_seconds <- times
  timings$reuse_mean_seconds <- mean(times)
  timings$reuse_median_seconds <- median(times)
  timings$reuse_max_estimate_difference <- max_est_diff
  timings$reuse_max_se_difference <- max_se_diff
} else if (a[4] != "primary") stop("mode must be primary or reuse")
rows <- lapply(seq_len(ncol(fit$params)), function(j) list(snp_index=j-1L,
  estimates=as.numeric(fit$params[,j]), ses=as.numeric(fit$SE_params[,j])))
out <- list(status="complete", mode=a[4], n=as.numeric(x$n), M=M, package=pkg,
            package_version=version, r_version=R.version.string, density="model", covariance_source=source,
            parameter_names=rownames(fit$params), timings=timings, rows=rows,
            script_elapsed_seconds=elapsed()-t_entry)
tmp <- paste0(a[3], ".tmp")
jsonlite::write_json(out, tmp, auto_unbox=TRUE, pretty=TRUE, digits=NA)
if (!file.rename(tmp, a[3])) stop("cannot atomically publish Stage 2 result")
