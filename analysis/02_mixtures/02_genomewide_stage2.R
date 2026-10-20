#!/usr/bin/env Rscript
# Genome-wide Stage 2 for the six mixture phenotypes
#
# What: maps every variant's 45 Stage 1 quantile effects onto the three mixture parameters
#   (gamma = logit P(lower component), beta_left, beta_right) by ordinary least squares with
#   the pooled Stage 1 covariance, under both weight designs, and adds a joint Hotelling test
#   of the three parameters. Writes the per-variant results to the working directory kept with the UK Biobank data
#   (they are read by the testosterone and antibody analyses) and exports genome-wide
#   summaries: genomic inflation, distance-clumped leads, Manhattan and QQ plot data.
#   The paper reports upper-component membership, -gamma; files keep gamma.
# Paper: Fig. 2; Fig. 3; Supplementary Fig. S9; Supplementary Tables S3, S4
# Input: Stage 1 quantile GWAS per chromosome (scan.stage1.tsv.gz, pooled covariance
#   scan.sigma0.tsv.gz and its metadata); weight matrices from 01_fit_mixtures.R
# Output: results/02_mixtures/stage2_mixtures/{fit_summary, leads_<fit>,
#   manhattan_points_<fit>_<design>, manhattan_bins_<fit>_<design>, qq_<fit>_<design>}.csv
# Runs: on a cluster with UK Biobank data. Chromosomes are independent and are run in
#   parallel on the cluster (Rscript 02_genomewide_stage2.R <chr> ...); the summary is written
#   once all 22 are done. About 3 minutes per phenotype, design and chromosome.

repo <- Sys.getenv("GMOQ_REPO", system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE))
source(file.path(repo, "R", "load.R"))
cfg <- read_config()
suppressPackageStartupMessages(library(gmoq))   # the paper used a development version of this package (Stage 2 commit 9fee67d)
stage2_commit <- cfg$package$stage2$commit
args <- commandArgs(trailingOnly = TRUE)
chromosomes <- if (length(args)) as.integer(args) else 1:22
n_threads <- as.integer(Sys.getenv("GMOQ_THREADS", "1"))
info_file <- function(design, fit, chr) sub("\\.tsv\\.gz$", ".info.csv", stage2_file(cfg, design, fit, chr))

## 1. Stage 2 for one phenotype, design and chromosome ------------------------------------
stage2_one <- function(f, design, chr, stage1_plain, stage1) {
  wobj <- readRDS(weight_file(cfg, design, f$fit, chr))
  W <- wobj$W; dir <- wobj$stage1_dir
  sigma0_file <- file.path(dir, "scan.sigma0.tsv.gz")
  nu <- as.numeric(read_json(file.path(dir, "scan.sigma0.meta.json"))$n_samples)
  beta_cols <- grep("^beta_tau_", names(stage1), value = TRUE)
  taus <- as.numeric(sub("^beta_tau_", "", beta_cols))

  # theta = A b with A = (W'W)^-1 W' for each variant's Stage 1 profile b; its covariance is
  # A Sigma0 A' / d, with d the residual genotype sum of squares (pooled covariance).
  t0 <- proc.time()[["elapsed"]]
  fit <- gmoq::param_gwas_from_file(
    stage1_file = stage1_plain, W = W, se_mode = "tau_cov", cov_source = "pooled",
    sigma0_file = sigma0_file, cov_meta_file = file.path(dir, "scan.cov.meta.json"),
    tau_cov_estimator = cfg$estimation$stage2_estimator, return_cov = "upper",
    param_cov_file = sub("\\.tsv\\.gz$", ".param_cov.gz", stage2_file(cfg, design, f$fit, chr)),
    n_threads = n_threads, q_method = "saddlepoint", require_density_mode = design)
  seconds <- proc.time()[["elapsed"]] - t0
  res <- as.data.table(if (is.list(fit) && !is.null(fit$results)) fit$results else fit)

  # Joint test of the three parameters: T = d theta' (A Sigma0 A')^-1 theta, referred to
  # Hotelling's T^2 with nu = the number of samples behind the pooled covariance.
  Wp <- unclass(W); attributes(Wp) <- attributes(Wp)[c("dim", "dimnames")]
  Sigma0 <- gw_read_sigma0(sigma0_file, taus)
  A <- solve(crossprod(Wp), t(Wp))
  V_unit <- A %*% Sigma0 %*% t(A)
  m <- match(res$snp_id, stage1$snp_id)
  if (anyNA(m)) stop("Stage 2 variants missing from Stage 1")
  B <- as.matrix(stage1[m, ..beta_cols]); d <- as.numeric(stage1$d_resid[m])
  theta <- B %*% t(A)
  joint_T <- d * rowSums((theta %*% solve(V_unit)) * theta)

  out <- data.table(snp_id = res$snp_id, chr = res$chr, bp = res$bp, effect_allele = res$effect_allele,
                    other_allele = res$other_allele, d_resid = d)
  for (p in colnames(Wp)) {
    out[[p]] <- res[[p]]; out[[paste0(p, "_se")]] <- res[[paste0(p, "_se")]]
    out[[paste0(p, "_p")]] <- res[[paste0(p, "_p")]]
    out[[paste0("mlog10p_", p)]] <- mlog10_p_hotelling((res[[p]] / res[[paste0(p, "_se")]])^2, 1, nu)
  }
  set(out, j = "joint_T", value = joint_T)
  set(out, j = "mlog10p_joint", value = mlog10_p_hotelling(joint_T, ncol(Wp), nu))
  fwrite(out, stage2_file(cfg, design, f$fit, chr), sep = "\t")
  fwrite(data.table(n_samples = nu, k = ncol(Wp), seconds = seconds,
                    density_weight_mode = gmoq::weight_density_mode(W, label = TRUE)),
         info_file(design, f$fit, chr))
  message(f$fit, " / ", design, " chr", chr, ": ", nrow(out), " variants, ", round(seconds, 1), " s")
}

for (chr in chromosomes) for (i in seq_len(nrow(MIXTURE_FITS))) {
  f <- MIXTURE_FITS[i]
  if (all(file.exists(vapply(WEIGHT_DESIGNS, function(d) info_file(d, f$fit, chr), "")))) next
  # The Stage 1 file is decompressed once and shared by the two designs.
  dir <- stage1_dir(cfg, f$cohort, chr)
  stage1_plain <- tempfile(fileext = ".tsv", tmpdir = private_dir(cfg, "tmp"))
  if (system2("gzip", c("-dc", file.path(dir, "scan.stage1.tsv.gz")), stdout = stage1_plain) != 0)
    stop("could not decompress ", dir)
  stage1 <- fread(stage1_plain, showProgress = FALSE)
  for (design in WEIGHT_DESIGNS) stage2_one(f, design, chr, stage1_plain, stage1)
  unlink(stage1_plain)
}

## 2. Genome-wide summaries ---------------------------------------------------------------
jobs <- CJ(fit = MIXTURE_FITS$fit, design = WEIGHT_DESIGNS, chr = 1:22)
if (!all(file.exists(jobs[, info_file(design, fit, chr), by = seq_len(nrow(jobs))]$V1))) { message("not all chromosomes are done; summary not written"); quit(save = "no") }

out_dir <- "02_mixtures/stage2_mixtures"
lead_p <- cfg$p07$lead_p                      # 5e-8
summary_rows <- list()
for (i in seq_len(nrow(MIXTURE_FITS))) {
  f <- MIXTURE_FITS[i]
  leads <- list()
  for (design in WEIGHT_DESIGNS) {
    info <- rbindlist(lapply(1:22, function(chr) fread(info_file(design, f$fit, chr))))
    gw <- read_stage2(cfg, design, f$fit)
    tests <- c(MIXTURE_PARAMS, "joint")
    mlog_cols <- paste0("mlog10p_", tests)
    nu <- info$n_samples[1]; k <- info$k[1]

    # Genomic inflation: median Wald statistic over the median of chi2(df).
    lambda <- vapply(tests, function(t) {
      stat <- if (t == "joint") gw$joint_T else (gw[[t]] / gw[[paste0(t, "_se")]])^2
      stats::median(stat, na.rm = TRUE) / qchisq(0.5, if (t == "joint") k else 1)
    }, numeric(1))
    n_gws <- vapply(mlog_cols, function(cc) sum(gw[[cc]] >= -log10(lead_p), na.rm = TRUE), numeric(1))
    row <- data.table(fit = f$fit, weight_design = design, trait = f$trait, model = f$model,
                      params = paste(MIXTURE_PARAMS, collapse = ","), n_chr = 22L, n_snps = nrow(gw),
                      n_samples = nu, nu = nu, k = k, density_weight_mode = info$density_weight_mode[1],
                      stage2_commit = stage2_commit, stage2_seconds_total = sum(info$seconds))
    for (t in tests) row[[paste0("lambda_gc_", t)]] <- lambda[[t]]
    for (t in seq_along(tests)) row[[paste0("n_gws_", tests[t])]] <- n_gws[[t]]
    summary_rows[[length(summary_rows) + 1L]] <- row

    # Leads: greedy 1 Mb distance clumping of each test at P < 5e-8.
    ld <- rbindlist(lapply(tests, function(t) {
      x <- as.data.table(distance_clump(gw, paste0("mlog10p_", t), lead_p, cfg$p07$lead_window_bp))
      if (nrow(x)) x[, lead_for_test := t]
      x
    }), fill = TRUE)
    if (nrow(ld)) {
      ld[, `:=`(fit = f$fit, weight_design = design)]
      setcolorder(ld, c("fit", "weight_design", "lead_for_test", "snp_id"))
      leads[[design]] <- ld
    }

    # Manhattan data (variants with any P < 1e-5, plus the maximum per 1 Mb bin) and QQ data.
    mt <- manhattan_thin(gw, mlog_cols, keep_p = cfg$p07$manhattan_keep_p, bin_bp = cfg$p07$manhattan_bin_bp)
    qq <- rbindlist(lapply(tests, function(t) {
      x <- as.data.table(qq_quantiles(gw[[paste0("mlog10p_", t)]])); x[, test := t]; x
    }))
    tag <- function(x) { x <- as.data.table(x); x[, weight_design := design]; setcolorder(x, "weight_design"); x }
    write_result(tag(mt$points), out_dir, paste0("manhattan_points_", f$fit, "_", design))
    write_result(tag(mt$bins), out_dir, paste0("manhattan_bins_", f$fit, "_", design))
    write_result(tag(qq), out_dir, paste0("qq_", f$fit, "_", design))
    rm(gw); invisible(gc())
  }
  if (length(leads)) write_result(rbindlist(leads, fill = TRUE), out_dir, paste0("leads_", f$fit))
}
write_result(rbindlist(summary_rows, fill = TRUE), out_dir, "fit_summary")
