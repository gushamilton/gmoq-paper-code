#!/usr/bin/env Rscript
# Interaction of each BMI lead variant with the BMI polygenic score
#
# What: for each of the 639 BMI lead variants, regresses BMI on the variant (dosage of the
#   BMI-increasing allele), the polygenic score that excludes the variant's chromosome, their
#   interaction, age, sex and ten genetic PCs, with heteroskedasticity-robust (HC0) and
#   model-based standard errors. Outcomes: log BMI (primary), BMI, rank-inverse-normal BMI and
#   the implied natural scale; the same fits are repeated in each random half for the
#   cross-fitted link test. For FTO (rs1421085) the variant's effect is also estimated within
#   each decile of the score. Finally each lead's half-sample Stage 1 profile is measured
#   against the score's own half-sample quantile profile. Everything written is per lead or per
#   decile; the summaries across leads are in 10_prs_summaries.R.
# Paper: Fig. 6a-c; Supplementary Fig. S14; Supplementary Table S6
# Input: the analysis sample and scores (07_prs_scores.R), the score's quantile profile
#   (08_prs_profile.py), lead dosages (plink2 .raw and .pvar per chromosome), the implied
#   natural-scale table, Stage 1 of the leads in each random half (LOCO none), and lead
#   profiles (results/01_data/lead_profiles)
# Output: results/05_bmi/prs_interaction/{per_lead_regressions,fto_prs_deciles,
#   prs_quantile_profile,outcome_summary}.csv
# Runs: on a cluster with UK Biobank data, in about 30 minutes

repo <- Sys.getenv("GMOQ_REPO", unset = NA)
if (is.na(repo)) repo <- system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE)
source(file.path(repo, "R", "load.R"))
source(file.path(repo, "analysis", "05_bmi", "bmi_common.R"))
cfg <- read_config()
P <- cfg$prs
u <- tau_grid(cfg)
PCS <- paste0("PC", 1:10)
private <- file.path(cfg$bp$private_root, "bmi_prs")
OUT <- "05_bmi/prs_interaction"

## 1. Inputs ------------------------------------------------------------------------------
A <- readRDS(file.path(private, "prep", "analysis.rds"))
n <- nrow(A); stopifnot(n == P$n)
profiles <- read_lead_profiles()
lead_set <- read_lead_set(profiles, P$lead_set, u, order_by_position = TRUE)   # log_bmi_prs_none
ref <- lead_set$meta[, .(snp_id, chr, bp, effect_allele, other_allele, d_resid)]

# Outcomes. The natural scale g(BMI) comes from a monotone spline of a tabulated transform,
# extended linearly beyond the tabulated range.
g_tab <- fread(P$natural_scale_table)
g_spline <- splinefun(g_tab$q_tau, g_tab$implied_scale_g, method = "monoH.FC")
q_lo <- min(g_tab$q_tau); q_hi <- max(g_tab$q_tau)
natural_scale <- function(q) {
  out <- g_spline(pmin(pmax(q, q_lo), q_hi))
  lo <- q < q_lo; hi <- q > q_hi
  out[lo] <- g_spline(q_lo) + g_spline(q_lo, deriv = 1) * (q[lo] - q_lo)
  out[hi] <- g_spline(q_hi) + g_spline(q_hi, deriv = 1) * (q[hi] - q_hi)
  out
}
Y <- cbind(raw = A$bmi, log = A$log_bmi, int = rank_inverse_normal(A$bmi), natural = natural_scale(A$bmi))
stopifnot(all(is.finite(Y)))
covariates <- as.matrix(A[, c("age", "sex", PCS), with = FALSE])
rows_of <- list(full = seq_len(n), h1 = which(A$half == 1), h2 = which(A$half == 2))

## 2. Per-lead interaction regressions ------------------------------------------------------
# y ~ snp + prs + snp:prs + age + sex + PC1-10. Kept per outcome: b (snp), g (interaction),
# b_prs (score main effect), and the (co)variances of b and g, HC0 (_h) and model-based (_m).
fit_interaction <- function(g, S, rows, outcomes) {
  X <- cbind(one = 1, snp = g[rows], prs = S[rows], int = g[rows] * S[rows], covariates[rows, , drop = FALSE])
  XtXi <- solve(crossprod(X))
  lapply(setNames(outcomes, outcomes), function(yc) {
    f <- ols_hc0(X, Y[rows, yc], XtXi)
    c(b = unname(f$b["snp"]), g = unname(f$b["int"]), b_prs = unname(f$b["prs"]),
      vb_h = f$V_hc0["snp", "snp"], vg_h = f$V_hc0["int", "int"], cbg_h = f$V_hc0["snp", "int"],
      vb_m = f$V_model["snp", "snp"], vg_m = f$V_model["int", "int"], cbg_m = f$V_model["snp", "int"])
  })
}
per_lead <- list(); fto <- list()
for (cc in sort(unique(ref$chr))) {
  pvar <- fread(file.path(P$lead_pvar_dir, sprintf("chr%02d_lead.pvar", cc)))
  setnames(pvar, 1:5, c("CHROM", "POS", "ID", "REF", "ALT"))
  dosage <- fread(file.path(P$lead_dosage_dir, sprintf("chr%02d_lead.raw", cc)), colClasses = list(character = c("FID", "IID")))
  row <- match(A$IID, dosage$IID); stopifnot(!anyNA(row))
  S <- A[[sprintf("PRS%02d", cc)]]
  for (i in which(ref$chr == cc)) {
    s <- ref$snp_id[i]; pv <- pvar[ID == s]
    # .raw columns are named <rsid>_<counted allele>(/<other allele>)
    col <- grep(paste0("^", s, "_[ACGT]+\\(/[ACGT]+\\)$"), names(dosage), value = TRUE)
    stopifnot(length(col) == 1L, nrow(pv) == 1L)
    counted <- sub("^[^_]+_([ACGT]+)\\(/([ACGT]+)\\)$", "\\1", col)
    other <- sub("^[^_]+_([ACGT]+)\\(/([ACGT]+)\\)$", "\\2", col)
    stopifnot(setequal(c(counted, other), c(pv$REF, pv$ALT)))
    x <- dosage[[col]][row]
    g <- if (ref$effect_allele[i] == counted) x else if (ref$effect_allele[i] == other) 2 - x else stop("allele mismatch ", s)
    out <- list(snp_id = s)
    full <- fit_interaction(g, S, rows_of$full, colnames(Y))
    for (yc in names(full)) for (k in names(full[[yc]])) out[[paste0(yc, "_", k)]] <- unname(full[[yc]][k])
    for (hh in c("h1", "h2")) {
      half <- fit_interaction(g, S, rows_of[[hh]], c("raw", "log"))
      for (yc in names(half)) for (k in c("b", "g", "vb_h", "vg_h", "cbg_h")) out[[paste0(hh, "_", yc, "_", k)]] <- unname(half[[yc]][k])
    }
    per_lead[[length(per_lead) + 1L]] <- as.data.table(out)
    # FTO: per-allele effect within each decile of the score (no interaction term)
    if (s == P$fto_snp) {
      decile <- pmin(pmax(findInterval(S, quantile(S, seq(0, 1, 0.1)), rightmost.closed = TRUE, all.inside = TRUE), 1L), 10L)
      for (dd in 1:10) for (yc in c("raw", "log")) {
        ix <- which(decile == dd)
        f <- ols_hc0(cbind(one = 1, snp = g[ix], covariates[ix, , drop = FALSE]), Y[ix, yc])
        fto[[length(fto) + 1L]] <- data.table(snp_id = s, scale = yc, prs_decile = dd, n = length(ix),
          mean_prs = mean(S[ix]), beta = unname(f$b["snp"]), se_hc0 = sqrt(f$V_hc0["snp", "snp"]),
          se_model = sqrt(f$V_model["snp", "snp"]))
      }
    }
  }
}
R <- rbindlist(per_lead)
stopifnot(identical(R$snp_id, ref$snp_id))

## 3. Half-sample profiles measured against the score's own profile -----------------------------
prs_profiles <- fread(file.path(private, "rif", "prs_profiles.csv"))
prs_profile_matrix <- function(trait, half) {
  x <- prs_profiles[scale == trait & sample == half][order(left_out_chr, tau)]
  matrix(x$b_prs, nrow = 22, byrow = TRUE)
}
# Stage 1 of the 639 leads in one random half (LOCO none), aligned to the reference allele
half_sample_leads <- function(arm) {
  dirs <- file.path(private, "stage1", arm, sprintf("chr%02d", 1:22))
  chrs <- sort(unique(ref$chr))
  X <- rbindlist(lapply(chrs, function(cc) {
    x <- fread(cmd = paste("zcat", shQuote(file.path(dirs[cc], "scan.stage1.tsv.gz"))), showProgress = FALSE)
    x[, snp_id := as.character(snp_id)][snp_id %in% ref$snp_id] }), fill = TRUE)
  beta_cols <- grep("^beta_tau_", names(X), value = TRUE)
  Sigma0 <- setNames(lapply(chrs, function(cc) gw_read_sigma0(file.path(dirs[cc], "scan.sigma0.tsv.gz"),
                                                               as.numeric(sub("^beta_tau_", "", beta_cols)))),
                     as.character(chrs))
  out <- align_to_reference(X[, .(snp_id, chr, bp, effect_allele, other_allele, d_resid)],
                            as.matrix(X[, beta_cols, with = FALSE]), ref)
  out$Sigma0 <- Sigma0
  out
}
for (trait in c("log_bmi", "bmi")) for (hh in c("h1", "h2")) {
  arm <- half_sample_leads(paste0(trait, "_", hh))
  tt <- prs_referenced_tilt(arm$B, arm$d, ref$chr, arm$Sigma0, prs_profile_matrix(trait, hh), u)
  for (k in names(tt)) set(R, j = paste("tilt", trait, hh, k, sep = "_"), value = tt[[k]])
}

## 4. Outputs ----------------------------------------------------------------------------
per_lead_regressions <- cbind(ref[, .(snp_id, chr, bp, effect_allele, other_allele)], R[, !"snp_id"])
outcome_summary <- data.table(scale = colnames(Y), n = n, mean = colMeans(Y), variance = apply(Y, 2, var))
write_result(per_lead_regressions, OUT, "per_lead_regressions")
write_result(rbindlist(fto), OUT, "fto_prs_deciles")
write_result(prs_profiles[sample == "full", .(scale, left_out_chr, tau, b_prs, se_prs, d_prs, n)], OUT, "prs_quantile_profile")
write_result(outcome_summary, OUT, "outcome_summary")
