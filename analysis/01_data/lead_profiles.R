#!/usr/bin/env Rscript
# Quantile-effect profiles of the lead variants, and the pooled covariance per chromosome
#
# What: pulls, for every lead set, the lead variants' Stage 1 rows from the genome-wide
#   quantile GWAS: the 45 quantile effects beta(u), their SEs and d_resid. Every lead is
#   expressed for the allele that increases location, delta_mu = sum_t c_t beta_t > 0.
#   BMI-scale sets that share the 639 BMI variants (log BMI, and the leads-only reruns used
#   for the polygenic-score analyses) use the allele that increases BMI in the main BMI
#   GWAS instead, so every BMI-scale table refers to the same allele. The pooled covariance
#   Sigma0 of each chromosome (and its mean over chromosomes) is exported with the
#   profiles: together with d_resid it gives the SE of any Stage 2 estimate,
#   Var(w'beta) = w' Sigma0 w / d_resid. Every downstream lead-level analysis starts here.
# Paper: Figs. 1-6; Supplementary Figs. S8-S14 (lead-level analyses)
# Input: lead sets (results/01_data/lead_variants/); Stage 1 quantile GWAS of each trait
#   (one file per chromosome); leads-only Stage 1 reruns on the 639 BMI leads without LOCO,
#   with LOCO as a covariate or subtracted, and by sex (analysis/00_gwas/03_stage1_quantile_gwas.sh)
# Output: results/01_data/lead_profiles/{lead_profiles,lead_profile_sets,rerun_reproduction,
#   sigma0_<source>}.csv
# Runs: on a cluster with the Stage 1 output, about 20 minutes (sets run 4 at a time)

.repo <- Sys.getenv("GMOQ_REPO", system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE))
source(file.path(.repo, "R", "load.R"))
cfg <- read_config()
u <- tau_grid(cfg)
cw <- trap_weights(u)
out <- "01_data/lead_profiles"

leads <- as.data.table(read_result("01_data/lead_variants", "leads"))
lead_sets <- as.data.table(read_result("01_data/lead_variants", "lead_sets"))

## 1. Which lead sets, from which Stage 1 GWAS
# The primary set of each trait, the sex-stratified testosterone sets and the weight
# delta_mu sets, each from its own trait's genome-wide Stage 1 run.
wanted <- lead_sets[primary == TRUE | rule == "sexstrat_B" | lead_set %in% c("weight_A_deltamu", "weight_B_deltamu")]
jobs <- lapply(seq_len(nrow(wanted)), function(i)
  list(lead_set = wanted$lead_set[i], trait = wanted$trait[i], source = wanted$trait[i],
       stage1 = stage1_chr_paths(cfg, wanted$trait[i])[stage1_exists == TRUE][order(chr)]))

# Leads-only Stage 1 reruns on the 639 BMI leads (same flags and sample as the genome-wide
# runs): raw BMI with no LOCO ("none", used for the polygenic-score work), LOCO as a
# covariate (as the genome-wide scan), LOCO subtracted, and within each sex; log BMI with
# no LOCO and with LOCO as a covariate.
bmi_ids <- leads[lead_set == "bmi_A_regenie"]$snp_id
leads_only_run <- function(dir, lead_set, trait) {
  d <- file.path(dir, sprintf("chr%02d", 1:22))
  ok <- file.exists(file.path(d, "scan.stage1.tsv.gz"))
  if (!any(ok)) return(NULL)
  list(lead_set = lead_set, trait = trait, source = lead_set, ids = bmi_ids,
       stage1 = data.table(chr = (1:22)[ok], stage1 = file.path(d, "scan.stage1.tsv.gz")[ok],
                           sigma0 = file.path(d, "scan.sigma0.tsv.gz")[ok]))
}
for (arm in c("none", "covariate", "subtract", "sex0", "sex1"))
  jobs[[length(jobs) + 1L]] <- leads_only_run(file.path(cfg$prs$stage1_leads_bmi, arm), paste0("bmi_prs_", arm), "bmi")
for (arm in c("none", "covariate"))
  jobs[[length(jobs) + 1L]] <- leads_only_run(file.path(cfg$prs$stage1_leads_log_bmi, paste0("log_bmi_", arm)),
                                              paste0("log_bmi_prs_", arm), "log_bmi")
jobs <- Filter(Negate(is.null), jobs)

## 2. Extract the lead rows from Stage 1
extract <- function(job) {
  ids <- job$ids %||% leads[lead_set == job$lead_set]$snp_id
  set_leads <- leads[lead_set == job$lead_set]
  files <- job$stage1
  # Read only the chromosomes that hold leads (all chromosomes when positions are unknown).
  chrs <- if (is.null(job$ids) && all(!is.na(set_leads$chr))) unique(set_leads$chr) else files$chr
  P <- rbindlist(lapply(intersect(files$chr, chrs), function(cc) {
    x <- read_stage1(files[chr == cc]$stage1, snps = ids)
    if (!nrow(x)) return(NULL)
    x[, snp_id := as.character(snp_id)]
  }), fill = TRUE)
  st <- stage1_taus(files$stage1[1])
  P <- P[, c("snp_id", "chr", "bp", "effect_allele", "other_allele", "d_resid", st$beta_cols, st$se_cols), with = FALSE]
  # A variant with two rows (multi-allelic record) is dropped.
  dup <- P[duplicated(snp_id)]$snp_id
  P <- P[!snp_id %in% dup]
  P[, delta_mu_raw := as.numeric(as.matrix(P[, st$beta_cols, with = FALSE]) %*% cw)]
  list(job = job, P = P, missing = setdiff(ids, P$snp_id), dup = unique(dup), st = st)
}
extracted <- parallel::mclapply(jobs, extract, mc.cores = 4, mc.preschedule = FALSE)
for (x in extracted) if (inherits(x, "try-error")) stop(x)

## 3. Align, add delta_mu and its pooled SE, and export Sigma0
# Reference allele for BMI-scale sets: the BMI-increasing allele in the main BMI GWAS.
bmi_main <- extracted[[which(vapply(extracted, function(x) x$job$lead_set, "") == "bmi_A_regenie")]]$P
bmi_ref <- bmi_main[, .(snp_id, ref_ea = effect_allele, ref_sign = sign(delta_mu_raw))]

profiles <- list(); set_rows <- list(); sigma0_written <- character()
for (x in extracted) {
  P <- copy(x$P); st <- x$st; job <- x$job
  bmi_scale <- job$lead_set == "log_bmi_A_regenie" || grepl("^(log_)?bmi_prs_", job$lead_set)
  if (bmi_scale) {
    m <- bmi_ref[match(P$snp_id, bmi_ref$snp_id)]
    same_ea <- P$effect_allele == m$ref_ea
    # flip unless this table's effect allele is the BMI-increasing allele
    flip <- !ifelse(same_ea, m$ref_sign >= 0, m$ref_sign < 0)
    own_disagree <- sum(sign(P$delta_mu_raw) != ifelse(flip, -1, 1), na.rm = TRUE)
  } else {
    flip <- P$delta_mu_raw < 0
    own_disagree <- 0L
  }
  flip[is.na(flip)] <- FALSE
  for (bc in st$beta_cols) set(P, which(flip), bc, -P[[bc]][flip])
  ea <- P$effect_allele; oa <- P$other_allele
  P[flip, `:=`(effect_allele = oa[flip], other_allele = ea[flip])]
  P[, aligned_flip := flip]

  # Pooled SE of delta_mu with each lead's own chromosome's Sigma0.
  Sigma0 <- lapply(split(job$stage1, by = "chr"), function(r) read_sigma0(r$sigma0, st$taus))
  var_factor <- vapply(Sigma0, function(M) as.numeric(t(cw) %*% M %*% cw), 0)
  P[, delta_mu := delta_mu_raw * ifelse(flip, -1, 1)]
  P[, delta_mu_se := sqrt(var_factor[as.character(chr)] / d_resid)]
  P[, delta_mu_z := delta_mu / delta_mu_se]
  P[, delta_mu_raw := NULL]
  P[, `:=`(lead_set = job$lead_set, trait = job$trait, stage1_source = job$source)]
  setcolorder(P, c("lead_set", "trait", "stage1_source", "snp_id", "chr", "bp", "effect_allele",
                   "other_allele", "aligned_flip", "d_resid", "delta_mu", "delta_mu_se", "delta_mu_z"))
  setorder(P, chr, bp)
  profiles[[job$lead_set]] <- P
  set_rows[[job$lead_set]] <- data.table(
    lead_set = job$lead_set, trait = job$trait, stage1_source = job$source,
    n_requested = length(job$ids %||% leads[lead_set == job$lead_set]$snp_id), n_found = nrow(P),
    n_missing = length(x$missing), missing_ids = paste(head(x$missing, 20), collapse = ";"),
    n_duplicate_dropped = length(x$dup), n_flipped = sum(flip),
    alignment = if (bmi_scale) "BMI-increasing allele (main BMI Stage 1)" else "own delta_mu > 0",
    n_own_sign_disagree = own_disagree, min_abs_z = min(abs(P$delta_mu_z)),
    median_d_resid = median(P$d_resid))

  # Sigma0 by chromosome and its mean, once per Stage 1 source: rows chr x tau, columns s_<tau>.
  if (!job$source %in% sigma0_written) {
    mean_Sigma0 <- Reduce(`+`, Sigma0) / length(Sigma0)
    as_rows <- function(M, label) data.table(chr = label, tau = st$taus, as.data.table(M))
    S <- rbindlist(c(lapply(names(Sigma0), function(k) as_rows(Sigma0[[k]], k)), list(as_rows(mean_Sigma0, "mean"))))
    setnames(S, c("chr", "tau", paste0("s_", sprintf("%.3f", st$taus))))
    write_result(S, out, paste0("sigma0_", job$source))
    sigma0_written <- c(sigma0_written, job$source)
  }
}

## 4. Leads-only reruns against the genome-wide profiles
# With LOCO as a covariate, the leads-only rerun must reproduce the genome-wide scan; the
# no-LOCO arm shows how much leaving out LOCO changes the profiles.
rerun_check <- function(a, b) {
  if (is.null(profiles[[a]]) || is.null(profiles[[b]])) return(NULL)
  bc <- grep("^beta_tau_", names(profiles[[a]]), value = TRUE)
  m <- merge(profiles[[a]][, c("snp_id", bc), with = FALSE], profiles[[b]][, c("snp_id", bc), with = FALSE], by = "snp_id")
  x <- as.matrix(m[, paste0(bc, ".x"), with = FALSE]); y <- as.matrix(m[, paste0(bc, ".y"), with = FALSE])
  data.table(lead_set = a, compared_to = b, n = nrow(m), max_abs_diff = max(abs(x - y)),
             profile_r = cor(as.numeric(x), as.numeric(y)))
}
reruns <- rbindlist(list(rerun_check("bmi_prs_covariate", "bmi_A_regenie"),
                         rerun_check("log_bmi_prs_covariate", "log_bmi_A_regenie"),
                         rerun_check("bmi_prs_none", "bmi_A_regenie"),
                         rerun_check("log_bmi_prs_none", "log_bmi_A_regenie")))
if (nrow(reruns)) write_result(reruns, out, "rerun_reproduction")

lead_profile_sets <- rbindlist(set_rows)
print(lead_profile_sets[, .(lead_set, n_requested, n_found, n_missing, n_flipped, n_own_sign_disagree)])
write_result(rbindlist(profiles, fill = TRUE), out, "lead_profiles")
write_result(lead_profile_sets, out, "lead_profile_sets")
