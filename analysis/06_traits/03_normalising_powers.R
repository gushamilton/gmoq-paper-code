#!/usr/bin/env Rscript
# Normalising and maximum-likelihood Box-Cox powers of each trait
#
# What: For each trait, the raw phenotype in exactly its Stage 1 analysis sample is summarised
#   by 03_boxcox_sample.py (aggregates only). The normalising power is the Box-Cox lambda that
#   maximises the normal probability-plot correlation of the transformed quantiles (193 levels
#   between the 2nd and 98th percentiles, lambda searched from -5 to 5 in steps of 0.01); the
#   textbook maximum-likelihood power and its 95% profile interval are reported beside it. The
#   pooled and per-lead flat powers from 02_flat_powers.R are joined on, so that one table
#   holds every power of a trait. The normalising power defines the Box-Cox scale used in
#   04_change_of_scale.R.
# Paper: Supplementary Table S5; Methods "Change of scale"
# Input: phenotype, covariate, keep and LOCO files of each trait's Stage 1 run (paths as listed
#   in results/01_data/sample_sizes/cohort_funnel.csv), the UK Biobank BGEN sample file;
#   results/06_traits/scale/{implied_scale,flat_lambda_summary}.csv
# Output: results/06_traits/scale/boxcox_lambda.csv
# Runs: on a cluster with UK Biobank data (one Python call per trait, a few minutes each)

repo <- Sys.getenv("GMOQ_REPO", system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE))
source(file.path(repo, "R", "load.R"))

cfg <- read_config()
lambda_grid <- round(seq(-500, 500) / 100, 2)
python <- cfg$package$stage1$python                     # Python with the pinned gmoq
helper <- file.path(repo, "analysis", "06_traits", "03_boxcox_sample.py")
work <- file.path(cfg$bp$private_root, "06_traits", "boxcox_sample")   # participant-level side stays here
dir.create(work, recursive = TRUE, showWarnings = FALSE)
cohorts <- as.data.table(read_result("01_data/sample_sizes", "cohort_funnel"))

survey <- cfg$cohorts$survey$traits
survey_traits <- names(survey)[vapply(survey, function(x) !identical(x$fig6, FALSE) &&
                                        !is.null(x[["n", exact = TRUE]]), logical(1))]
traits <- c("bmi", "weight", survey_traits)

## 1. Aggregates of each raw trait in its Stage 1 sample
sample_summary <- function(trait_name) {
  r <- cohorts[trait == trait_name]
  out <- file.path(work, paste0(trait_name, "_boxcox.json"))
  log_file <- sub("json$", "log", out)
  args <- c(helper, "--trait", trait_name, "--pheno", r$pheno_file, "--pheno-col", r$pheno_col,
            "--covar", r$covar_file, "--covar-cols", shQuote(r$covariates), "--keep", r$keep_file,
            "--loco", if (is.na(r$loco_file_chr1) || !nzchar(r$loco_file_chr1)) "none" else r$loco_file_chr1,
            # any one chromosome's BGEN supplies the sample order
            "--bgen", file.path(cfg$ukb$bgen_dir, "data.chr22.bgen"), "--sample", cfg$ukb$bgen_sample,
            "--out", out)
  status <- system2(python, args, stdout = log_file, stderr = log_file)
  if (status != 0L || !file.exists(out)) stop("Box-Cox sample summary failed for ", trait_name)
  jsonlite::read_json(out, simplifyVector = TRUE)
}
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
summaries <- parallel::mclapply(traits, sample_summary, mc.cores = min(n_cores, 4L), mc.preschedule = FALSE)
names(summaries) <- traits

## 2. Flat powers from 02_flat_powers.R
implied <- as.data.table(read_result("06_traits/scale", "implied_scale"))
flat_summary <- as.data.table(read_result("06_traits/scale", "flat_lambda_summary"))
# Pooled flat power: the Box-Cox curve that best matches the implied natural scale g(u).
pooled_flat_power <- function(trait_name) {
  x <- implied[trait == trait_name]
  fit_r <- vapply(lambda_grid, function(l) if (any(x$q_raw <= 0)) NA_real_ else
    suppressWarnings(cor(boxcox(x$q_raw, l), x$g)), 0)
  if (all(is.na(fit_r))) NA_real_ else lambda_grid[which.max(fit_r)]
}

## 3. Normalising and maximum-likelihood powers
powers <- rbindlist(lapply(traits, function(tr) {
  s <- summaries[[tr]]
  p <- s$p_grid; qy <- s$q_raw_type7
  lambda_main <- ppcor_main <- NA_real_
  if (all(qy > 0)) {
    ppcor <- vapply(lambda_grid, function(l) cor(boxcox(qy, l), qnorm(p)), 0)
    lambda_main <- lambda_grid[which.max(ppcor)]; ppcor_main <- max(ppcor)
  }
  lambda_ml <- ml_lo <- ml_hi <- NA_real_
  if (!is.null(s$boxcox_loglik)) {
    ll <- s$boxcox_loglik
    lambda_ml <- lambda_grid[which.max(ll)]
    ci <- range(lambda_grid[ll >= max(ll) - qchisq(0.95, 1) / 2])   # 95% profile-likelihood interval
    ml_lo <- ci[1]; ml_hi <- ci[2]
  }
  data.table(trait = tr, analysed_scale = if (grepl("^log_", tr)) "log" else "raw", n = s$n,
             n_nonpositive = s$n_nonpositive, skew_raw = s$skew_raw,
             skew_log = if (is.null(s$skew_log)) NA_real_ else s$skew_log,
             lambda_main = lambda_main, ppcor_main = ppcor_main, ppcor_raw = cor(qy, qnorm(p)),
             ppcor_log = if (all(qy > 0)) cor(log(qy), qnorm(p)) else NA_real_,
             lambda_main_at_grid_edge = !is.na(lambda_main) && abs(abs(lambda_main) - 5) < 1e-9,
             lambda_ml = lambda_ml, lambda_ml_lo95 = ml_lo, lambda_ml_hi95 = ml_hi,
             lambda_flat_pooled = pooled_flat_power(tr),
             lambda_v_median_strong = flat_summary[trait == tr]$lambda_v_median_strong)
}))
print(powers)

## 4. Write
write_result(powers, "06_traits/scale", "boxcox_lambda")
