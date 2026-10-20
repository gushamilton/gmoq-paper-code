#!/usr/bin/env Rscript
# Genome-wide Stage 2 for the BMI models, with heritabilities and genetic correlations
#
# What: applies four BMI models (additive, multiplicative about zero, additive +
#   multiplicative, and log-scale location + variance) to the Stage 1 quantile GWAS of BMI at
#   every variant, giving a GWAS of each model parameter. It then estimates the SNP
#   heritability of each parameter, and the genetic correlations among them and with the
#   standard REGENIE BMI GWAS, by LD score regression. The jointly estimated additive and
#   multiplicative coefficients are nearly collinear, so their mutual genetic correlation
#   depends on the origin of the multiplicative term and is not interpreted.
# Paper: Methods (BMI models, genome-wide)
# Input: Stage 1 quantile GWAS of log BMI (or BMI) with its pooled covariance, one file per
#   chromosome; the standard REGENIE GWAS of the same trait; HapMap3 SNP list and 1000
#   Genomes European LD scores; model weights (results/05_bmi/model_scores)
# Output: results/05_bmi/genomewide_models/{fit_summary,lambda_gc,weight_matrices,h2,
#   rg_matrix}.csv, manhattan_bins_<model>.csv, qq_<model>.csv
# Runs: on a cluster with UK Biobank data (chromosomes can be run in parallel; about an hour
#   per chromosome, then LDSC). The paper used a development version of the
#   G-MoQ R package (Stage 2 commit 9fee67d); gmoq 0.7.0 has the same estimators.

repo <- Sys.getenv("GMOQ_REPO", unset = NA)
if (is.na(repo)) repo <- system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE)
source(file.path(repo, "R", "load.R"))
source(file.path(repo, "analysis", "05_bmi", "bmi_common.R"))
suppressPackageStartupMessages(library(gmoq))   # G-MoQ Stage 2
cfg <- read_config()
COHORT <- cfg$headline$bmi_scale                         # "log_bmi" (paper) or "bmi"
SCALE <- c(log_bmi = "log", bmi = "raw")[[COHORT]]
MODELS <- unlist(cfg$p12b$models)                        # additive, log_additive, both, logscale_vqtl
LEAD_P <- cfg$p12b$lead_p; BIN_BP <- cfg$p12b$manhattan_bin_bp
ldsc <- cfg$p08                                          # LDSC module, HapMap3 list, LD score directory
u <- tau_grid(cfg)
OUT <- "05_bmi/genomewide_models"
scratch <- file.path(cfg$bp$private_root, "05_bmi_genomewide", COHORT)   # per-variant files stay here
for (s in c("stage1", "raw", "munged", "h2", "rg")) dir.create(file.path(scratch, s), recursive = TRUE, showWarnings = FALSE)

## 1. Model weights -----------------------------------------------------------------------
# Exactly the matrices used for the lead-variant scores, on this scale.
definitions <- as.data.table(read_result("05_bmi/model_scores", "model_definitions"))[scale == SCALE]
weight_table <- as.data.table(read_result("05_bmi/model_scores", "weight_matrices"))[trait == "bmi" & scale == SCALE]
W_of <- lapply(setNames(MODELS, MODELS), function(m) {
  cols <- strsplit(definitions[model == m]$columns, ",")[[1]]
  x <- weight_table[model == m][order(as.numeric(tau))]
  W <- as.matrix(x[, cols, with = FALSE]); storage.mode(W) <- "double"
  rownames(W) <- sprintf("%.3f", u)
  stopifnot(all(is.finite(W)), check_rank(W)$full)          # every parameter needs its own signature
  W
})

## 2. Stage 2 at every variant, chromosome by chromosome --------------------------------------
# OLS theta = A beta with A = (W'W)^-1 W', and its pooled covariance A Sigma0 A' / d. The
# package fits; a direct computation checks it. Tests are Hotelling T with nu = the number of
# samples behind Sigma0 (1 df per parameter, k df jointly).
fits <- stage2_seconds <- setNames(vector("list", length(MODELS)), MODELS)
for (chr in 1:22) {
  stage1_dir <- file.path(cfg$cohorts[[COHORT]]$root, sprintf("chr%02d", chr))
  stage1_file <- file.path(scratch, "stage1", sprintf("chr%02d.tsv", chr))
  stopifnot(system2("gzip", c("-dc", file.path(stage1_dir, "scan.stage1.tsv.gz")), stdout = stage1_file) == 0)
  s1 <- fread(stage1_file, showProgress = FALSE)
  beta_cols <- grep("^beta_tau_", names(s1), value = TRUE)
  stopifnot(isTRUE(all.equal(as.numeric(sub("^beta_tau_", "", beta_cols)), round(u, 3), tolerance = 5e-7)))
  sigma0_file <- file.path(stage1_dir, "scan.sigma0.tsv.gz")
  Sigma0 <- gw_read_sigma0(sigma0_file, u)
  nu <- as.numeric(jsonlite::read_json(file.path(stage1_dir, "scan.sigma0.meta.json"), simplifyVector = TRUE)$n_samples)
  stopifnot(nu == cfg$cohorts[[COHORT]]$n)
  for (m in MODELS) {
    W <- W_of[[m]]
    t0 <- proc.time()[["elapsed"]]
    fit <- gmoq::param_gwas_from_file(
      stage1_file = stage1_file, W = W, se_mode = "tau_cov", cov_source = "pooled",
      sigma0_file = sigma0_file, cov_meta_file = file.path(stage1_dir, "scan.cov.meta.json"),
      tau_cov_estimator = "ols", return_cov = "none", n_threads = 1L, q_method = "saddlepoint",
      require_density_mode = "model")
    seconds <- proc.time()[["elapsed"]] - t0
    res <- as.data.table(if (is.list(fit) && !is.null(fit$results)) fit$results else fit)
    A <- solve(crossprod(W), t(W)); V_unit <- A %*% Sigma0 %*% t(A)
    row <- match(res$snp_id, s1$snp_id); stopifnot(!anyNA(row))
    B <- as.matrix(s1[row, ..beta_cols]); d <- as.numeric(s1$d_resid[row])
    theta <- B %*% t(A); params <- colnames(W)
    for (j in seq_along(params)) stopifnot(
      max(abs(res[[params[j]]] - theta[, j])) <= 1e-6 * max(1e-12, max(abs(theta[, j]))),
      max(abs(res[[paste0(params[j], "_se")]] / sqrt(V_unit[j, j] / d) - 1)) <= 1e-6)
    out <- data.table(snp_id = res$snp_id, chr = res$chr, bp = res$bp,
                      effect_allele = res$effect_allele, other_allele = res$other_allele)
    for (p in params) {
      out[[p]] <- res[[p]]; out[[paste0(p, "_se")]] <- res[[paste0(p, "_se")]]
      out[[paste0("mlog10p_", p)]] <- mlog10_p_hotelling((res[[p]] / res[[paste0(p, "_se")]])^2, 1, nu)
    }
    if (length(params) > 1) {
      out[, joint_T := d * rowSums((theta %*% solve(V_unit)) * theta)]
      out[, mlog10p_joint := mlog10_p_hotelling(joint_T, length(params), nu)]
    }
    out[, n_samples := nu]
    fits[[m]] <- rbind(fits[[m]], out)
    stage2_seconds[[m]] <- c(stage2_seconds[[m]], seconds)
  }
  unlink(stage1_file)
  message("chr", chr, " done")
}

## 3. Genome-wide summaries per model ------------------------------------------------------
# munge_sumstats + ldsc --h2 for one GWAS (signed BETA, P, N), returning its munged file.
sh <- function(cmd) {                                       # login shell so that `module` exists
  out <- system2("bash", c("-lc", shQuote(cmd)), stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(out, "status")) && attr(out, "status") != 0) stop("command failed: ", cmd, "\n", paste(tail(out, 20), collapse = "\n"))
  out
}
munge <- function(id, dt) {
  raw <- file.path(scratch, "raw", paste0(id, ".tsv"))
  fwrite(dt[is.finite(BETA) & is.finite(SE) & SE > 0 & is.finite(P)], raw, sep = "\t")
  sh(sprintf("module load %s && munge_sumstats.py --sumstats %s --out %s --merge-alleles %s --snp SNP --a1 A1 --a2 A2 --p P --N-col N --signed-sumstats BETA,0 --ignore SE",
             ldsc$ldsc_module, raw, file.path(scratch, "munged", id), ldsc$ldsc_merge_alleles))
  unlink(raw)
  file.path(scratch, "munged", paste0(id, ".sumstats.gz"))
}
fit_summary <- lambda_rows <- weight_rows <- list(); ldsc_files <- list()
for (m in MODELS) {
  gw <- fits[[m]]; params <- colnames(W_of[[m]]); k <- length(params)
  tests <- c(params, if (k > 1) "joint")
  p_cols <- paste0("mlog10p_", tests)
  # lambda_GC = median chi-square / its null median
  lambda <- vapply(tests, function(t) {
    stat <- if (t == "joint") gw$joint_T else (gw[[t]] / gw[[paste0(t, "_se")]])^2
    median(stat, na.rm = TRUE) / qchisq(0.5, if (t == "joint") k else 1)
  }, 0)
  n_gws <- vapply(p_cols, function(cc) sum(gw[[cc]] >= -log10(LEAD_P), na.rm = TRUE), 0)
  lambda_rows[[m]] <- data.table(model = m, test = tests, lambda_gc = lambda, n_snps = nrow(gw), n_gws = n_gws,
                                 df = ifelse(tests == "joint", k, 1))
  row <- data.table(model = m, scale = SCALE, cohort = COHORT, params = paste(params, collapse = ","),
                    n_chr = 22L, n_snps = nrow(gw), n_samples = gw$n_samples[1], k = k,
                    stage2_commit = cfg$package$stage2$commit, stage2_seconds_total = sum(stage2_seconds[[m]]))
  for (t in seq_along(tests)) { row[[paste0("lambda_gc_", tests[t])]] <- lambda[[t]]; row[[paste0("n_gws_", tests[t])]] <- n_gws[[t]] }
  fit_summary[[m]] <- row
  for (j in seq_along(params)) weight_rows[[length(weight_rows) + 1L]] <- data.table(model = m, scale = SCALE,
    tau = u, param = params[j], W = as.numeric(W_of[[m]][, j]))
  # Manhattan plot bins (maximum -log10 p per 1 Mb) and QQ quantiles: aggregate only
  manhattan <- manhattan_thin(gw, p_cols, keep_p = 1e-300, bin_bp = BIN_BP)
  write_result(cbind(data.table(model = m), manhattan$bins), OUT, paste0("manhattan_bins_", m))
  write_result(rbindlist(lapply(tests, function(t) cbind(data.table(model = m, test = t),
    as.data.table(qq_quantiles(gw[[paste0("mlog10p_", t)]]))))), OUT, paste0("qq_", m))
  for (p in params) {
    id <- paste0(m, "__", p)
    ldsc_files[[id]] <- list(model = m, param = p, file = munge(id, data.table(SNP = gw$snp_id, CHR = gw$chr, BP = gw$bp,
      A1 = gw$effect_allele, A2 = gw$other_allele, N = gw$n_samples, P = pmax(10^(-gw[[paste0("mlog10p_", p)]]), 1e-300),
      BETA = gw[[p]], SE = gw[[paste0(p, "_se")]])))
  }
}

## 4. Standard REGENIE GWAS of the same trait ---------------------------------------------
regenie_paths <- regenie_files(cfg, COHORT, cfg$p12b$regenie_role)[order(chr)]$path
regenie <- rbindlist(lapply(regenie_paths, function(f) fread(f, select = c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1",
                                                                         "N", "TEST", "BETA", "SE", "LOG10P"))[TEST == "ADD"]))
regenie_id <- paste0("regenie__", COHORT)
ldsc_files <- c(setNames(list(list(model = "regenie", param = "BETA", file = munge(regenie_id,
  regenie[, .(SNP = ID, CHR = CHROM, BP = GENPOS, A1 = ALLELE1, A2 = ALLELE0, N, P = pmax(10^(-LOG10P), 1e-300), BETA, SE)]))),
  regenie_id), ldsc_files)
rm(regenie)

## 5. LD score regression ------------------------------------------------------------------
pair_of_numbers <- function(lines, pattern) {
  x <- grep(pattern, lines, value = TRUE)[1]
  if (is.na(x)) return(c(NA_real_, NA_real_))
  v <- regmatches(sub("^[^:]*:", "", x), gregexpr("-?[0-9.]+(e-?[0-9]+)?", sub("^[^:]*:", "", x)))[[1]]
  c(suppressWarnings(as.numeric(v[1])), if (length(v) > 1) suppressWarnings(as.numeric(v[2])) else NA_real_)
}
parse_h2_log <- function(lg) {
  h2 <- pair_of_numbers(lg, "^Total Observed scale h2"); lgc <- pair_of_numbers(lg, "^Lambda GC")
  mc <- pair_of_numbers(lg, "^Mean Chi"); ic <- pair_of_numbers(lg, "^Intercept"); rt <- pair_of_numbers(lg, "^Ratio:")
  n_snps <- suppressWarnings(as.integer(sub(".*After merging with regression SNP LD, ([0-9]+) SNPs remain.*", "\\1",
              grep("After merging with regression SNP LD", lg, value = TRUE)[1])))
  data.table(h2 = h2[1], h2_se = h2[2], ldsc_lambda_gc = lgc[1], mean_chi2 = mc[1], intercept = ic[1],
             intercept_se = ic[2], ratio = rt[1], ratio_se = rt[2],
             ratio_note = if (any(grepl("^Ratio <", lg))) trimws(grep("^Ratio <", lg, value = TRUE)[1]) else NA_character_,
             n_snps_ldsc = n_snps)
}
parse_rg_log <- function(lg) {                              # the summary table at the end of an rg log
  i <- grep("^\\s*p1\\s+p2\\s+rg", lg)[1]
  header <- strsplit(trimws(lg[i]), "\\s+")[[1]]
  rows <- list(); j <- i + 1L
  while (j <= length(lg) && nzchar(trimws(lg[j]))) {
    v <- strsplit(trimws(lg[j]), "\\s+")[[1]]
    if (length(v) == length(header)) rows[[length(rows) + 1L]] <- as.data.table(setNames(as.list(v), header))
    j <- j + 1L
  }
  rbindlist(rows)
}
ids <- names(ldsc_files)
h2 <- rbindlist(lapply(ids, function(id) {
  prefix <- file.path(scratch, "h2", id)
  sh(sprintf("module load %s && ldsc.py --h2 %s --ref-ld-chr %s --w-ld-chr %s --out %s",
             ldsc$ldsc_module, ldsc_files[[id]]$file, ldsc$ldsc_ld_dir, ldsc$ldsc_ld_dir, prefix))
  cbind(data.table(trait = id, model = ldsc_files[[id]]$model, param = ldsc_files[[id]]$param),
        parse_h2_log(readLines(paste0(prefix, ".log"))))
}))
# Genetic correlations for all pairs: trait i against traits i+1..n in one ldsc call, run from
# the munged directory with file basenames so that the summary table has one row per pair.
rg <- rbindlist(lapply(seq_len(length(ids) - 1L), function(i) {
  prefix <- file.path(scratch, "rg", ids[i]); others <- ids[(i + 1L):length(ids)]
  sh(sprintf("module load %s && cd %s && ldsc.py --rg %s --ref-ld-chr %s --w-ld-chr %s --out %s", ldsc$ldsc_module,
             shQuote(file.path(scratch, "munged")),
             paste(basename(vapply(c(ids[i], others), function(o) ldsc_files[[o]]$file, "")), collapse = ","),
             ldsc$ldsc_ld_dir, ldsc$ldsc_ld_dir, prefix))
  x <- parse_rg_log(readLines(paste0(prefix, ".log")))
  stopifnot(nrow(x) == length(others))
  id_of <- function(p) sub("\\.sumstats\\.gz$", "", basename(p))
  num <- function(v) suppressWarnings(as.numeric(v))
  x[, .(trait1 = id_of(p1), trait2 = id_of(p2), rg = num(rg), rg_se = num(se), z = num(z), p = num(p),
        gcov_int = num(gcov_int), gcov_int_se = num(gcov_int_se))]
}))
rg[, `:=`(h2_1 = h2$h2[match(trait1, h2$trait)], h2_2 = h2$h2[match(trait2, h2$trait)])]
setcolorder(rg, c("trait1", "trait2", "rg", "rg_se", "p", "h2_1", "h2_2", "gcov_int"))

write_result(rbindlist(fit_summary, fill = TRUE), OUT, "fit_summary")
write_result(rbindlist(lambda_rows), OUT, "lambda_gc")
write_result(rbindlist(weight_rows), OUT, "weight_matrices")
write_result(h2, OUT, "h2")
write_result(rg, OUT, "rg_matrix")
