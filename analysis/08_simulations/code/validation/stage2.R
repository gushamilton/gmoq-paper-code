#!/usr/bin/env Rscript
# Stage 2 for one simulated variant
#
# What: reads one variant's Stage 1 summary (45 quantile effects, their pooled or jackknife
#   covariance, the sample quantiles and the null-model fit) written by run_examples.py,
#   builds the model's weight matrix W (normal: location and dispersion; two-component
#   mixture: membership and the two component locations) with gmoq, and fits Stage 2
#   with gmoq::param_gwas (OLS, the chosen covariance). It also returns the exact
#   population projection target passed in from truth.py. It is a thin adapter around the
#   package and does not reimplement the estimator.
# Paper: Supplementary Figs. S3-S6; Supplementary Table S2 (studies 3-7)
# Input: --input one-variant Stage 1 JSON; --cov-source pooled|jackknife; the package
#   checkout, commit and R library of the environment made by bootstrap.sh
# Output: --output one-variant Stage 2 JSON (estimates, SEs, confidence intervals, P values,
#   target)
# Runs: on a cluster, called once per variant and covariance by the drivers, about a second

PINNED_PACKAGE_COMMIT <- "9c47e7ba56370e708810fff8c2c9ad3b723cbccf"
PINNED_PACKAGE_VERSION <- "0.4.0"
PINNED_STAGE1_VERSION <- "0.1.0"

parse_cli <- function(args) {
  if (length(args) == 0L || length(args) %% 2L != 0L) {
    stop(paste(
      "Usage: Rscript --vanilla validation/stage2.R",
      "--input INPUT.json --output OUTPUT.json",
      "--cov-source pooled|jackknife --package-root DIR",
      "--package-commit SHA --package-version VERSION --r-lib DIR",
      "--environment-manifest FILE"
    ), call. = FALSE)
  }
  keys <- args[seq(1L, length(args), by = 2L)]
  vals <- args[seq(2L, length(args), by = 2L)]
  if (any(!grepl("^--", keys))) stop("Stage 2 arguments must use --name value", call. = FALSE)
  out <- as.list(vals)
  names(out) <- sub("^--", "", keys)
  out
}

arg_required <- function(args, name) {
  value <- args[[name]]
  if (is.null(value) || !nzchar(value)) stop("Missing --", name, call. = FALSE)
  value
}

as_numeric_vector <- function(x, name) {
  value <- suppressWarnings(as.numeric(unlist(x, use.names = FALSE)))
  if (!length(value) || any(!is.finite(value))) stop(name, " is not finite", call. = FALSE)
  value
}

as_named_list <- function(x, names_x) {
  values <- as.numeric(x)
  if (length(values) != length(names_x)) stop("unexpected vector length", call. = FALSE)
  as.list(stats::setNames(values, names_x))
}

row_major_upper <- function(S) {
  S <- as.matrix(S)
  if (nrow(S) != ncol(S)) stop("covariance matrix is not square", call. = FALSE)
  idx <- which(upper.tri(S, diag = TRUE))
  ord <- order(row(S)[idx], col(S)[idx])
  as.numeric(S[idx[ord]])
}

check_source_pin <- function(package_root, expected_sha) {
  if (!dir.exists(package_root)) stop("Pinned package root does not exist: ", package_root, call. = FALSE)
  sha <- tryCatch(
    system2("git", c("-C", package_root, "rev-parse", "HEAD"), stdout = TRUE, stderr = TRUE),
    error = function(e) stop("Unable to verify pinned source with git", call. = FALSE)
  )
  sha <- trimws(as.character(sha))
  if (!length(sha) || !any(sha == expected_sha)) {
    stop("Refusing Stage 2 source: package root is not at pinned commit ", expected_sha,
         " (git reported: ", paste(sha, collapse = " | "), ")", call. = FALSE)
  }
  tracked <- tryCatch(
    system2("git", c("-C", package_root, "status", "--porcelain=v1", "--untracked-files=no"),
            stdout = TRUE, stderr = TRUE),
    error = function(e) stop("Unable to inspect pinned source with git", call. = FALSE)
  )
  if (!is.null(attr(tracked, "status")) && attr(tracked, "status") != 0L) {
    stop("git status failed for pinned package source", call. = FALSE)
  }
  tracked <- trimws(as.character(tracked))
  if (length(tracked) && any(nzchar(tracked))) {
    stop("Refusing Stage 2 source: tracked package files are modified: ",
         paste(tracked, collapse = " | "), call. = FALSE)
  }
  expected_sha
}

manifest_string <- function(value, label) {
  if (is.null(value) || length(value) != 1L || is.na(value) ||
      !nzchar(as.character(value))) stop("Environment manifest is missing ", label, call. = FALSE)
  as.character(value)
}

manifest_path <- function(value, label) {
  normalizePath(manifest_string(value, label), mustWork = TRUE)
}

check_environment_manifest <- function(path, r_lib, package_root) {
  manifest <- tryCatch(
    jsonlite::fromJSON(path, simplifyVector = FALSE),
    error = function(e) stop("Cannot parse environment manifest: ", conditionMessage(e), call. = FALSE)
  )
  if (!is.list(manifest) || !is.list(manifest$package) || !is.list(manifest$stage2) ||
      !is.list(manifest$r) || !is.list(manifest$r$gmoq)) {
    stop("Environment manifest does not match bootstrap schema", call. = FALSE)
  }
  manifest_root <- manifest_path(manifest$package$root, "package.root")
  if (!identical(manifest_root, package_root)) stop("--package-root differs from manifest package.root", call. = FALSE)
  if (!identical(manifest_string(manifest$package$commit, "package.commit"), PINNED_PACKAGE_COMMIT) ||
      !identical(manifest_string(manifest$package$head_verified, "package.head_verified"), PINNED_PACKAGE_COMMIT)) {
    stop("Environment manifest does not record the pinned package commit", call. = FALSE)
  }
  if (!isTRUE(manifest$package$tracked_clean_before_build) ||
      !isTRUE(manifest$package$tracked_clean_after_build)) {
    stop("Environment manifest does not certify clean tracked package source", call. = FALSE)
  }
  manifest_rlib <- manifest_path(manifest$stage2$isolated_library, "stage2.isolated_library")
  recorded_rlib <- manifest_path(manifest$r$gmoq$isolated_library,
                                 "r.gmoq.isolated_library")
  if (!identical(manifest_rlib, recorded_rlib) || !identical(manifest_rlib, r_lib)) {
    stop("--r-lib is not the exact private Rlib recorded in the environment manifest", call. = FALSE)
  }
  installed <- manifest_path(manifest$r$gmoq$installed_path,
                             "r.gmoq.installed_path")
  if (!startsWith(installed, paste0(manifest_rlib, .Platform$file.sep))) {
    stop("Manifest gmoq path is outside its private Rlib", call. = FALSE)
  }
  if (!identical(manifest_string(manifest$r$gmoq$version,
                                 "r.gmoq.version"), PINNED_PACKAGE_VERSION)) {
    stop("Environment manifest records the wrong gmoq version", call. = FALSE)
  }
  description_sha <- tolower(manifest_string(manifest$stage2$description_sha256,
                                              "stage2.description_sha256"))
  if (!grepl("^[[:xdigit:]]{64}$", description_sha)) {
    stop("Environment manifest has an invalid Stage 2 DESCRIPTION SHA", call. = FALSE)
  }
  list(rlib = manifest_rlib, installed_path = installed, description_sha = description_sha)
}

sha256_file <- function(path) {
  if (!nzchar(Sys.which("sha256sum"))) stop("sha256sum is required to verify installed DESCRIPTION", call. = FALSE)
  output <- tryCatch(
    system2("sha256sum", c("--", path), stdout = TRUE, stderr = TRUE),
    error = function(e) stop("Unable to hash installed DESCRIPTION with sha256sum", call. = FALSE)
  )
  if (!length(output) || (!is.null(attr(output, "status")) && attr(output, "status") != 0L)) {
    stop("sha256sum failed for installed DESCRIPTION", call. = FALSE)
  }
  hash <- sub("[[:space:]].*$", "", trimws(as.character(output[[1L]])))
  if (!grepl("^[[:xdigit:]]{64}$", hash)) stop("sha256sum returned an invalid hash", call. = FALSE)
  tolower(hash)
}

check_weight_design <- function(W, label) {
  rank <- qr(W)$rank
  condition <- suppressWarnings(kappa(W, exact = TRUE))
  if (rank < ncol(W) || !is.finite(condition) || condition > 1e12) {
    stop(label, " is rank deficient or ill-conditioned: rank=", rank, "/", ncol(W),
         ", condition=", signif(condition, 6), call. = FALSE)
  }
  list(rank = as.integer(rank), n_columns = ncol(W), condition_number = condition,
       condition_number_limit = 1e12)
}

make_weights <- function(family, taus, q_tau, baseline) {
  if (identical(family, "normal")) {
    gmoq::make_weights_vqtl(
      taus = taus,
      q_tau = q_tau,
      mu = as.numeric(baseline$mu),
      sd = as.numeric(baseline$sd)
    )
  } else if (identical(family, "mixture")) {
    gmoq::make_weights_normal_mixture(
      taus = taus,
      q_tau = q_tau,
      p1 = as.numeric(baseline$p1),
      mu1 = as.numeric(baseline$mu1),
      sd1 = as.numeric(baseline$sd1),
      mu2 = as.numeric(baseline$mu2),
      sd2 = as.numeric(baseline$sd2),
      include_membership = TRUE
    )
  } else {
    stop("Unsupported family: ", family, call. = FALSE)
  }
}

run_one <- function(input, cov_source, installed_stage2_path) {
  family <- as.character(input$family)
  # run_examples.py writes a flat, one-SNP aggregate hand-off, not a nested
  # input$stage1 object. Require the schema rather than accepting empty vectors.
  required_fields <- c("taus", "q_tau", "baseline_fit", "d_resid", "population_truth")
  if (!all(required_fields %in% names(input))) {
    stop("Stage 1 hand-off is missing required top-level fields", call. = FALSE)
  }
  stage1_in <- input
  taus <- as_numeric_vector(stage1_in$taus, "stage1 taus")
  q_sample <- as_numeric_vector(stage1_in$q_tau, "stage1 q_tau")
  if (length(taus) != length(q_sample)) stop("tau and q_tau lengths differ", call. = FALSE)
  beta_name <- if (identical(cov_source, "pooled")) "beta_pooled" else "beta_jackknife"
  se_name <- if (identical(cov_source, "pooled")) "se_tau_pooled" else "se_tau_jackknife"
  beta <- as_numeric_vector(stage1_in[[beta_name]], beta_name)
  se_tau <- as_numeric_vector(stage1_in[[se_name]], se_name)
  if (length(beta) != length(taus) || length(se_tau) != length(taus)) {
    stop("Stage 1 tau arrays do not have the configured length", call. = FALSE)
  }
  d_resid <- as.numeric(stage1_in$d_resid)
  if (length(d_resid) != 1L || !is.finite(d_resid) || d_resid <= 0) {
    stop("Stage 1 d_resid must be one positive finite value", call. = FALSE)
  }

  # Production W is always made by the pinned R package.  W_hat uses the
  # same-sample nuisance fit and the sample Stage 1 q_tau, as required by the
  # validation contract. W_pop is only the fixed-target diagnostic counterpart.
  W_hat <- make_weights(family, taus, q_sample, stage1_in$baseline_fit)
  truth <- input$population_truth
  q_pop <- as_numeric_vector(truth$q_tau, "population q_tau")
  baseline_pop <- truth$baseline
  W_pop <- make_weights(family, taus, q_pop, baseline_pop)
  if (nrow(W_hat) != length(taus) || nrow(W_pop) != length(taus)) {
    stop("package weight builder returned the wrong number of rows", call. = FALSE)
  }
  weight_diagnostics_hat <- check_weight_design(W_hat, "same-sample fitted W")
  weight_diagnostics_pop <- check_weight_design(W_pop, "population W")
  W_independent <- as.matrix(truth$W_independent)
  if (!all(dim(W_independent) == dim(W_pop)) || any(!is.finite(W_independent))) {
    stop("independent truth weights have the wrong shape or are non-finite", call. = FALSE)
  }
  weight_error <- max(abs(W_pop - W_independent))

  if (identical(cov_source, "pooled")) {
    sigma0 <- as.matrix(stage1_in$sigma0)
    if (!all(dim(sigma0) == c(length(taus), length(taus))) || any(!is.finite(sigma0))) {
      stop("pooled Sigma_0 has the wrong shape or is non-finite", call. = FALSE)
    }
    sigma <- 0.5 * (sigma0 + t(sigma0)) / d_resid
    covariance_upper <- row_major_upper(sigma)
  } else {
    covariance_upper <- as_numeric_vector(
      stage1_in$cov_upper_jackknife, "jackknife covariance"
    )
    expected <- length(taus) * (length(taus) + 1L) / 2L
    if (length(covariance_upper) != expected) {
      stop("packed jackknife covariance has the wrong length", call. = FALSE)
    }
  }
  if (any(!is.finite(covariance_upper))) stop("covariance is non-finite", call. = FALSE)

  # param_gwas consumes the same cov_data layout as the package file reader:
  # packed row-major upper triangles and byte offsets.  One SNP therefore has
  # one payload beginning at byte offset zero. No hand-written estimator is used.
  stage1 <- list(
    taus = taus,
    q_tau = q_sample,
    Q_slope = matrix(beta, nrow = length(taus), ncol = 1L,
                     dimnames = list(paste0("tau", taus), "sim_snp_1")),
    SE_tau = matrix(se_tau, nrow = length(taus), ncol = 1L,
                    dimnames = list(paste0("tau", taus), "sim_snp_1"))
  )
  cov_data <- list(cov_vec = covariance_upper, offsets = 0)
  is_pooled <- identical(cov_source, "pooled")
  result <- gmoq::param_gwas(
    stage1 = stage1,
    transform = "custom_W",
    transform_args = list(W = W_hat),
    se_mode = "tau_cov",
    cov_data = cov_data,
    tau_cov_estimator = "ols",
    n_threads = 1L,
    return_cov = "upper",
    q_method = "saddlepoint",
    nu = if (is_pooled) as.numeric(input$n) else NULL,
    nu_source = if (is_pooled) "pooled_sigma0_n_samples" else NULL,
    n_blocks = if (is_pooled) NULL else as.integer(input$n_blocks),
    n_eff = if (is_pooled) NULL else d_resid,
    unknown_nu = "error"
  )

  parameter_names <- colnames(W_hat)
  if (is.null(parameter_names)) parameter_names <- paste0("parameter_", seq_len(ncol(W_hat)))
  theta <- as.numeric(result$params[, 1L])
  se <- as.numeric(result$SE_params[, 1L])
  p <- if (is.null(result$P_params)) NULL else as.numeric(result$P_params[, 1L])
  reference_distribution <- as.character(result$param_ref_dist)
  if (length(reference_distribution) != 1L || is.na(reference_distribution) ||
      !nzchar(reference_distribution)) {
    stop("param_gwas returned a missing or non-scalar reference distribution", call. = FALSE)
  }
  df_matrix <- if (is.null(result$param_ref_df)) NULL else as.matrix(result$param_ref_df)
  df <- if (is.null(df_matrix)) NULL else as.numeric(df_matrix[, 1L])
  if (length(theta) != length(parameter_names) || length(se) != length(theta) ||
      any(!is.finite(theta)) || any(!is.finite(se)) || any(se < 0)) {
    stop("param_gwas returned invalid parameter estimates or standard errors", call. = FALSE)
  }
  critical <- if (identical(reference_distribution, "t")) {
    if (is.null(df_matrix) || !all(dim(df_matrix) == dim(result$params)) ||
        length(df) != length(theta) || any(!is.finite(df)) || any(df <= 0)) {
      stop("param_gwas returned missing or invalid t-reference degrees of freedom", call. = FALSE)
    }
    stats::qt(0.975, df = df)
  } else if (identical(reference_distribution, "normal")) {
    rep(stats::qnorm(0.975), length(theta))
  } else {
    stop("Cannot construct confidence intervals for reference distribution: ",
         reference_distribution, call. = FALSE)
  }
  ci_lower <- theta - critical * se
  ci_upper <- theta + critical * se
  b_pop <- as_numeric_vector(truth$b_tau, "population b_tau")
  if (length(b_pop) != length(taus)) stop("population b_tau length differs", call. = FALSE)
  theta_pop <- as_numeric_vector(truth$theta, "population theta")
  if (length(theta_pop) != ncol(W_hat)) stop("population theta length differs from W", call. = FALSE)
  A_hat <- result$A
  theta_random_weight <- as.numeric(A_hat %*% b_pop)
  A_pop <- solve(crossprod(W_pop), t(W_pop))
  theta_pop_recomputed <- as.numeric(A_pop %*% b_pop)
  if (length(theta_pop_recomputed) != length(theta_pop)) stop("population target dimension mismatch", call. = FALSE)

  estimate <- list(
    target_kind = "fixed_population_theta",
    parameter_names = parameter_names,
    estimate = as_named_list(theta, parameter_names),
    se = as_named_list(se, parameter_names),
    p_value = if (is.null(p)) NULL else as_named_list(p, parameter_names),
    df = if (is.null(df)) NULL else as_named_list(df, parameter_names),
    ci_lower = as_named_list(ci_lower, parameter_names),
    ci_upper = as_named_list(ci_upper, parameter_names),
    Q = as.numeric(result$Q[1L]),
    Q_p_value = if (is.null(result$Q_pval)) NULL else as.numeric(result$Q_pval[1L]),
    Q_df = as.numeric(result$df),
    reference_distribution = reference_distribution,
    reference_source = as.character(result$param_ref_source),
    covariance_source = cov_source,
    tau_cov_estimator = "ols"
  )
  list(
    status = "ok",
    family = family,
    covariance_source = cov_source,
    inference_role = if (is_pooled) "primary" else "jackknife_sensitivity",
    estimator = "OLS",
    n = as.integer(input$n),
    n_blocks = if (is_pooled) NULL else as.integer(input$n_blocks),
    stage1_density_method = as.character(input$density_method),
    stage1_density_convention = as.character(input$density_convention),
    model_weight_convention = as.character(input$stage2_weight_convention),
    same_sample_fitted_baseline = input$baseline_fit,
    weights = list(
      fitted_W = W_hat,
      population_W = W_pop,
      fitted_diagnostics = weight_diagnostics_hat,
      population_diagnostics = weight_diagnostics_pop,
      independent_population_W_max_abs_error = weight_error
    ),
    target = list(
      kind = "fixed_population_theta",
      theta_pop = as_named_list(theta_pop, parameter_names),
      theta_pop_from_package_W = as_named_list(theta_pop_recomputed, parameter_names),
      generating = truth$generating,
      diagnostics = truth$diagnostics
    ),
    random_weight_diagnostic = list(
      kind = "A_W_hat_times_b_pop_not_fixed_target",
      value = as_named_list(theta_random_weight, parameter_names)
    ),
    parameter_covariance_upper = result$param_cov_upper,
    estimate = estimate,
    stage1_package = input$package,
    package = list(
      name = "gmoq",
      version = as.character(utils::packageVersion("gmoq")),
      location = installed_stage2_path
    )
  )
}

args <- parse_cli(commandArgs(trailingOnly = TRUE))
input_path <- arg_required(args, "input")
output_path <- arg_required(args, "output")
cov_source <- arg_required(args, "cov-source")
if (!(cov_source %in% c("pooled", "jackknife"))) stop("--cov-source must be pooled or jackknife", call. = FALSE)
package_root <- normalizePath(arg_required(args, "package-root"), mustWork = TRUE)
expected_sha <- arg_required(args, "package-commit")
expected_version <- arg_required(args, "package-version")
r_lib <- normalizePath(arg_required(args, "r-lib"), mustWork = TRUE)
environment_manifest <- normalizePath(arg_required(args, "environment-manifest"), mustWork = TRUE)
if (!identical(expected_sha, PINNED_PACKAGE_COMMIT)) stop("wrong compact package commit argument", call. = FALSE)
if (!identical(expected_version, PINNED_PACKAGE_VERSION)) stop("wrong compact package version argument", call. = FALSE)
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("validation/stage2.R requires jsonlite; install it in the R library", call. = FALSE)
}
environment <- check_environment_manifest(environment_manifest, r_lib, package_root)
source_sha <- check_source_pin(package_root, PINNED_PACKAGE_COMMIT)
.libPaths(unique(c(environment$rlib, .libPaths())))
stage2_path <- find.package("gmoq", lib.loc = environment$rlib, quiet = TRUE)
if (!length(stage2_path)) stop("gmoq is not installed in the manifest private Rlib", call. = FALSE)
stage2_path <- normalizePath(stage2_path[[1L]], mustWork = TRUE)
if (!startsWith(stage2_path, paste0(environment$rlib, .Platform$file.sep)) ||
    !identical(stage2_path, environment$installed_path)) {
  stop("installed gmoq resolved outside the manifest installation path", call. = FALSE)
}
library(gmoq, lib.loc = environment$rlib)
if (!identical(as.character(utils::packageVersion("gmoq")), PINNED_PACKAGE_VERSION)) {
  stop("installed gmoq is not version ", PINNED_PACKAGE_VERSION, call. = FALSE)
}
description_path <- normalizePath(file.path(stage2_path, "DESCRIPTION"), mustWork = TRUE)
if (!identical(sha256_file(description_path), environment$description_sha)) {
  stop("installed gmoq DESCRIPTION SHA does not match environment manifest", call. = FALSE)
}
input <- jsonlite::fromJSON(input_path, simplifyVector = TRUE)
n_samples <- suppressWarnings(as.integer(input$n))
n_blocks <- suppressWarnings(as.integer(input$n_blocks))
if (length(n_samples) != 1L || is.na(n_samples) ||
    length(n_blocks) != 1L || is.na(n_blocks) ||
    !(n_samples %in% c(2000L, 10000L, 50000L, 400000L)) || n_samples < n_blocks) {
  stop("compact Stage 2 adapter supports N=2000, N=10000, N=50000 or N=400000 with N >= n_blocks", call. = FALSE)
}
if (is.null(input$density_method) || !identical(as.character(input$density_method), "silverman_gaussian")) {
  stop("compact Stage 2 adapter requires the silverman_gaussian Stage 1 density estimate", call. = FALSE)
}
if (is.null(input$package$package_commit) || !identical(as.character(input$package$package_commit), source_sha)) {
  stop("input Stage 1 provenance does not carry the verified pinned package commit", call. = FALSE)
}
if (is.null(input$package$stage1_distribution_version) ||
    !identical(as.character(input$package$stage1_distribution_version), PINNED_STAGE1_VERSION)) {
  stop("input Stage 1 provenance does not carry the pinned distribution version ", PINNED_STAGE1_VERSION, call. = FALSE)
}
if (file.exists(output_path)) stop("refusing to overwrite Stage 2 output: ", output_path, call. = FALSE)
output <- run_one(input, cov_source, stage2_path)
output$source_check <- list(
  package_root = package_root,
  package_commit = source_sha,
  package_version = PINNED_PACKAGE_VERSION,
  stage1_distribution_version = PINNED_STAGE1_VERSION,
  stage2_source_checked = TRUE,
  compiled_stage1_required = TRUE,
  numpy_fallback_allowed = FALSE
)
# Preserve numerical inference precision; jsonlite's default rounds to 4 digits.
jsonlite::write_json(output, output_path, auto_unbox = TRUE, pretty = TRUE,
                     na = "null", digits = NA)
