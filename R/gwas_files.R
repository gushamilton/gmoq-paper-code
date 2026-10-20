# Locating and reading the per-chromosome GWAS files (Stage 1 and REGENIE)
#
# What: functions used by every script that reads the genome-wide Stage 1 output: which
#   traits have a Stage 1 GWAS, where each chromosome's Stage 1 table, pooled covariance and
#   log are, and readers for those files and for REGENIE step 2 output. Loaded by R/load.R.
# Input: config.yaml (cohorts: Stage 1 roots; regenie_gwas: REGENIE step 2 files)
# Runs: wherever the GWAS files are (a cluster with the UK Biobank-derived files)
#
# Stage 1 output per chromosome (gmoq-stage1 --out <prefix>):
#   <prefix>.stage1.tsv.gz        one row per variant: snp_id, chr, bp, alleles, d_resid,
#                                 beta_tau_<u> and se_tau_<u> at the 45 quantile levels
#   <prefix>.sigma0.tsv.gz        pooled 45 x 45 covariance Sigma0 (the SE of any linear
#                                 combination a'beta is sqrt(a' Sigma0 a / d_resid))
#   <prefix>.sigma0.meta.json     N, covariates, LOCO mode, kernel
#   stage1.log                    sample counts at each filter
#   rif_nuisance.json             (some runs) quantiles, kernel density and mean of the trait

suppressPackageStartupMessages(library(data.table))
`%||%` <- function(a, b) if (is.null(a)) b else a

# A cohort's N in config.yaml. The YAML key `n` can be read as the boolean FALSE by
# YAML 1.1 parsers, so both spellings are accepted.
cfg_n <- function(x) {
  v <- x[["n", exact = TRUE]]
  if (is.null(v)) v <- x[["FALSE", exact = TRUE]]
  if (is.null(v)) NA_real_ else as.numeric(v)
}

# Traits with a Stage 1 GWAS: the main cohorts, then the trait-survey traits except the
# vitamin D run that did not complete (config note "no vitamin D").
stage1_traits <- function(cfg) {
  main <- setdiff(names(cfg$cohorts), "survey")
  sv <- cfg$cohorts$survey$traits
  surv <- names(sv)[!vapply(sv, function(x) identical(x$note, "no vitamin D"), logical(1))]
  data.table(
    trait = c(main, surv),
    group = c(rep("main", length(main)), rep("survey", length(surv))),
    root = c(vapply(main, function(t) cfg$cohorts[[t]]$root, ""),
             file.path(cfg$cohorts$survey$root, surv)),
    n_config = c(vapply(main, function(t) cfg_n(cfg$cohorts[[t]]), 0),
                 vapply(surv, function(t) cfg_n(sv[[t]]), 0)),
    fig6 = c(rep(NA, length(main)),
             vapply(surv, function(t) !identical(sv[[t]]$fig6, FALSE), logical(1))))
}

# One row per chromosome: where the Stage 1 table, pooled covariance and logs are.
# Most runs are laid out as <root>/chrNN/scan.*. The main cohorts' runs record each
# chromosome's files in <run>/progress.json instead, because some chromosomes were taken
# from an earlier identical run (same package, LOCO mode, kernel, N and covariance).
stage1_chr_paths <- function(cfg, trait) {
  tt <- stage1_traits(cfg)
  root <- tt$root[match(trait, tt$trait)]
  if (is.na(root)) stop("unknown trait: ", trait)
  progress <- file.path(dirname(dirname(root)), "progress.json")
  rows <- NULL
  if (file.exists(progress)) {
    pr <- jsonlite::read_json(progress, simplifyVector = FALSE)
    # progress.json writes testosterone_shrink_0.2; config writes testosterone_shrink_0_2
    ch <- Filter(function(e) identical(gsub(".", "_", e$trait, fixed = TRUE), trait), pr$chromosomes)
    if (length(ch)) rows <- rbindlist(lapply(ch, function(e) data.table(
      chr = as.integer(e$chr), stage1 = e$stage1, cov_prefix = e$cov_prefix,
      source_kind = e$source_kind %||% NA_character_,
      package_commit = e$package_commit %||% NA_character_,
      progress_n = as.numeric(e$n %||% NA))))
  }
  if (is.null(rows)) {
    d <- file.path(root, sprintf("chr%02d", 1:22))
    rows <- data.table(chr = 1:22, stage1 = file.path(d, "scan.stage1.tsv.gz"),
                       cov_prefix = file.path(d, "scan"), source_kind = "root_chr_dir",
                       package_commit = NA_character_, progress_n = NA_real_)
  }
  setorder(rows, chr)
  rows[, `:=`(trait = trait, dir = dirname(stage1),
              sigma0 = paste0(cov_prefix, ".sigma0.tsv.gz"),
              sigma0_meta = paste0(cov_prefix, ".sigma0.meta.json"))]
  rows[, stage1_log := file.path(dir, "stage1.log")]
  rows[, rif_nuisance := ifelse(file.exists(file.path(dir, "rif_nuisance.json")),
                                file.path(dir, "rif_nuisance.json"), NA_character_)]
  # A zero-byte table is an interrupted run, not a result.
  rows[, stage1_exists := file.exists(stage1) & file.size(stage1) > 0]
  rows[]
}

# Sample counts and input files from a Stage 1 log (never returns raw log lines, some of
# which list example sample identifiers).
parse_stage1_log <- function(path) {
  out <- list(log_found = file.exists(path))
  if (!out$log_found) return(out)
  x <- readLines(path, warn = FALSE)
  grab <- function(re, num = TRUE) {
    m <- Filter(length, regmatches(x, regexec(re, x)))
    if (!length(m)) return(if (num) NA_real_ else NA_character_)
    if (num) as.numeric(m[[1]][2]) else m[[1]][2]
  }
  out$keep_file <- grab("Loading sample keep list from (\\S+)", num = FALSE)
  out$n_keep <- grab("Keeping ([0-9]+) samples")
  out$pheno_file <- grab("Loading phenotypes from (\\S+)", num = FALSE)
  out$covar_file <- grab("Loading covariates from (\\S+)", num = FALSE)
  out$n_valid_pheno <- grab("Loaded ([0-9]+) samples with valid phenotypes")
  out$loco_file <- grab("Loading LOCO predictions for \\S+ from (\\S+)", num = FALSE)
  out$n_absent_loco <- grab("\\] ([0-9]+) phenotype samples are absent from LOCO")
  out$n_bgen_matched <- grab("Matched ([0-9]+)/[0-9]+ samples to BGEN")
  out$n_before_bgen <- grab("Matched [0-9]+/([0-9]+) samples to BGEN")
  out
}

# Column names and quantile levels of a Stage 1 table.
stage1_taus <- function(path) {
  hdr <- names(fread(cmd = sprintf("zcat %s | head -1", shQuote(path)), header = TRUE))
  bc <- grep("^beta_tau_", hdr, value = TRUE)
  list(beta_cols = bc, se_cols = sub("^beta_", "se_", bc), taus = as.numeric(sub("^beta_tau_", "", bc)))
}

# Stage 1 rows, optionally only the variants in `snps` (filtered with awk before reading).
read_stage1 <- function(path, cols = NULL, snps = NULL) {
  if (is.null(snps))
    return(fread(cmd = sprintf("zcat %s", shQuote(path)), select = cols, showProgress = FALSE))
  ids <- tempfile(fileext = ".txt"); writeLines(unique(snps), ids); on.exit(unlink(ids))
  x <- fread(cmd = sprintf("zcat %s | awk 'NR==FNR{k[$1]=1;next} FNR==1 || ($1 in k)' %s -",
                           shQuote(path), shQuote(ids)), header = TRUE, showProgress = FALSE)
  if (!is.null(cols)) x <- x[, intersect(cols, names(x)), with = FALSE]
  x
}

# Pooled Sigma0 as a symmetric 45 x 45 matrix in the order of `taus`.
read_sigma0 <- function(path, taus = NULL) {
  s0 <- fread(path, showProgress = FALSE)
  sc <- grep("^s_", names(s0), value = TRUE)
  if (is.null(taus)) taus <- as.numeric(s0$tau)
  ri <- match(round(taus, 6), round(as.numeric(s0$tau), 6))
  ci <- match(round(taus, 6), round(as.numeric(sub("^s_", "", sc)), 6))
  if (anyNA(ri) || anyNA(ci)) stop("Sigma0 quantile levels do not match: ", path)
  M <- as.matrix(s0[, sc, with = FALSE])[ri, ci, drop = FALSE]
  M <- (M + t(M)) / 2
  dimnames(M) <- list(sprintf("%.3f", taus), sprintf("%.3f", taus))
  M
}

# REGENIE step 2 files of a trait for a role ("comparator": the additive GWAS that the
# Stage 1 mean is compared with; "leads": the GWAS that is clumped for lead variants), from
# config.yaml regenie_gwas. A path with {chr} is one file per chromosome (chr = 01..22);
# otherwise one genome-wide file.
regenie_files <- function(cfg, trait, role) {
  p <- cfg$regenie_gwas[[trait]][[role]]
  if (is.null(p)) return(data.table(chr = integer(), path = character(), layout = character()))
  if (grepl("{chr}", p, fixed = TRUE)) {
    data.table(chr = 1:22, path = vapply(sprintf("%02d", 1:22), function(cc) gsub("{chr}", cc, p, fixed = TRUE), ""),
               layout = "per_chr")
  } else data.table(chr = NA_integer_, path = p, layout = "combined")
}

# REGENIE step 2 additive rows, optionally one chromosome of a genome-wide file.
read_regenie <- function(path, chr = NULL) {
  cols <- c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "N", "TEST", "BETA", "SE", "LOG10P")
  cmd <- if (grepl("\\.gz$", path)) sprintf("zcat %s", shQuote(path)) else sprintf("cat %s", shQuote(path))
  if (!is.null(chr)) cmd <- sprintf("%s | awk 'NR==1 || $1==%d'", cmd, as.integer(chr))
  x <- fread(cmd = cmd, select = cols, showProgress = FALSE)
  if ("TEST" %in% names(x)) x <- x[TEST == "ADD"]
  x[, CHROM := as.integer(CHROM)]
  x[]
}

# Sign that aligns a REGENIE effect (ALLELE1 is the effect allele) to the Stage 1
# effect_allele: +1 same alleles, -1 swapped, NA otherwise (no strand flips).
allele_sign <- function(ea, oa, a1, a0) {
  s <- rep(NA_real_, length(ea))
  s[ea == a1 & oa == a0] <- 1
  s[ea == a0 & oa == a1] <- -1
  s
}
