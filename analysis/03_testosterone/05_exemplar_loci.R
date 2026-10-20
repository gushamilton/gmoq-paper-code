#!/usr/bin/env Rscript
# Exemplar testosterone variants: one per kind of effect
#
# What: picks five lead variants of the unshrunk testosterone analysis that illustrate the
#   kinds of effect a two-component model can separate, by a rule fixed before selection.
#   The pool is the distance-clumped (1 Mb) leads of the joint test at P < 5e-8; z = estimate/SE.
#     female_specific    largest |z_left|  with |z_right| < 1.96
#     male_specific      largest |z_right| with |z_left| < 1.96
#     shared_concordant  largest min(|z_left|, |z_right|), same sign, both >= 5.45
#     discordant         largest min(|z_left|, |z_right|), opposite signs, both >= 5.45
#     membership         largest |z_gamma| with |z_left| < 1.96 and |z_right| < 1.96
#   Ties go to the larger joint -log10 P. The rule is applied under each weight design; the
#   variants chosen under the reported design for unshrunk testosterone (hybrid) are then shown
#   in every arm and design: their Stage 1 profile, its fit W theta and each parameter's part
#   W_j theta_j, the estimates with the sex-stratified REGENIE effects, and the GWAS of sex at
#   each exemplar aligned to the same allele. gamma is logit P(lower = female component); the
#   sex GWAS codes male = 1, so it is compared with -gamma (upper-component membership).
# Paper: Fig. 2d,e
# Input: per-variant Stage 2 results and weight matrices (02_mixtures); Stage 1 profiles and
#   pooled covariance; sex-stratified REGENIE GWAS; the sex GWAS (config.yaml, section p08)
# Output: results/03_testosterone/recovery/{exemplar_selected,exemplar_profiles,exemplar_estimates}.csv;
#   results/03_testosterone/exemplars/exemplar_sex_reference.csv
# Runs: on a cluster with UK Biobank-derived Stage 1 output, in about 20 minutes

repo <- Sys.getenv("GMOQ_REPO", system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE))
source(file.path(repo, "R", "load.R"))
cfg <- read_config()
p8 <- cfg$p08
ARMS <- p8$arms
PRIMARY <- mixture_design(cfg, "testosterone")                    # hybrid
DESIGNS <- c(PRIMARY, setdiff(WEIGHT_DESIGNS, PRIMARY))

# Sex-stratified REGENIE GWAS (as in 03_compare_with_sex_stratified_gwas.R).
read_combined <- function(arm, sex) {
  x <- fread(file.path(p8$regenie_combined, sprintf("testosterone_%s_%s_all.regenie", sex, arm)),
             select = c("ID", "ALLELE0", "ALLELE1", "N", "BETA", "SE", "LOG10P"), showProgress = FALSE)
  data.table(snp_id = x$ID, a1 = x$ALLELE1, a2 = x$ALLELE0, n = x$N, beta = x$BETA, se = x$SE, mlog10p = x$LOG10P)
}
read_ldsc_format <- function(path) {
  x <- fread(path, select = c("SNP", "A1", "A2", "N", "P", "BETA", "SE"), showProgress = FALSE)
  data.table(snp_id = x$SNP, a1 = x$A1, a2 = x$A2, n = x$N, beta = x$BETA, se = x$SE, mlog10p = -log10(pmax(x$P, 1e-320)))
}
REF <- list()
for (a in ARMS) for (s in c("female", "male")) {
  REF[[paste("combined", a, s)]] <- read_combined(a, s)
  REF[[paste("paper_matched", a, s)]] <- read_ldsc_format(file.path(p8$regenie_paper_raw, sprintf("regenie__%s_%s.tsv", a, s)))
}

selected <- list(); profiles <- list(); estimates <- list(); exemplar_ids <- NULL
for (design in DESIGNS) {
  G <- setNames(lapply(ARMS, function(a) read_stage2(cfg, design, a)), ARMS)

  ## 1. Apply the selection rule to unshrunk testosterone --------------------------------
  g <- copy(G[["testosterone"]])
  g[, `:=`(zL = beta_left / beta_left_se, zR = beta_right / beta_right_se, zG = gamma / gamma_se)]
  pool <- as.data.table(distance_clump(g, "mlog10p_joint", cfg$p07$lead_p, cfg$p07$lead_window_bp))
  pick <- function(sub, score) if (!nrow(sub)) NULL else sub[order(-score, -mlog10p_joint)][1]
  sel <- list(
    female_specific = pick(pool[abs(zR) < 1.96], pool[abs(zR) < 1.96, abs(zL)]),
    male_specific = pick(pool[abs(zL) < 1.96], pool[abs(zL) < 1.96, abs(zR)]),
    shared_concordant = { s <- pool[sign(zL) == sign(zR) & pmin(abs(zL), abs(zR)) >= 5.45]; pick(s, s[, pmin(abs(zL), abs(zR))]) },
    discordant = { s <- pool[sign(zL) != sign(zR) & pmin(abs(zL), abs(zR)) >= 5.45]; pick(s, s[, pmin(abs(zL), abs(zR))]) },
    membership = { s <- pool[abs(zL) < 1.96 & abs(zR) < 1.96]; pick(s, s[, abs(zG)]) })
  ex <- rbindlist(lapply(names(sel), function(n) if (is.null(sel[[n]])) data.table(category = n) else
    cbind(data.table(category = n), sel[[n]][, .(snp_id, chr, bp)])), fill = TRUE)
  ex[, pool_size := nrow(pool)]
  selected[[design]] <- cbind(data.table(weight_design = design), ex)
  if (design == PRIMARY) exemplar_ids <- ex[!is.na(snp_id)]

  ## 2. Profiles and estimates of the reported-design exemplars, every arm ---------------
  for (a in ARMS) for (i in seq_len(nrow(exemplar_ids))) {
    e <- exemplar_ids[i]
    wobj <- readRDS(weight_file(cfg, design, a, e$chr)); W <- wobj$W
    s1 <- read_stage1(wobj$stage1_dir)[snp_id == e$snp_id]
    beta_cols <- grep("^beta_tau_", names(s1), value = TRUE)
    taus <- as.numeric(sub("^beta_tau_", "", beta_cols))
    Sigma0 <- gw_read_sigma0(file.path(wobj$stage1_dir, "scan.sigma0.tsv.gz"), taus)
    beta <- as.numeric(s1[1, ..beta_cols]); d <- s1$d_resid[1]
    gg <- G[[a]][snp_id == e$snp_id]
    theta <- as.numeric(gg[1, colnames(W), with = FALSE])
    pr <- data.table(weight_design = design, arm = a, category = e$category, snp_id = e$snp_id, tau = taus,
                     beta_stage1 = beta, se_stage1 = sqrt(diag(Sigma0) / d), fitted_total = as.numeric(W %*% theta))
    for (j in seq_along(theta)) pr[[paste0("fitted_", colnames(W)[j])]] <- W[, j] * theta[j]
    profiles[[length(profiles) + 1L]] <- pr
    er <- gg[, .(weight_design = design, arm = a, category = e$category, snp_id, chr, bp, effect_allele, other_allele,
                 gamma, gamma_se, beta_left, beta_left_se, beta_right, beta_right_se, mlog10p_joint)]
    for (s in c("female", "male")) for (rn in c("combined", "paper_matched")) {
      h <- harmonise_with_reference(G[[a]][snp_id == e$snp_id], "beta_left", 1, REF[[paste(rn, a, s)]])
      er[[paste0("regenie_", rn, "_", s, "_beta")]] <- if (nrow(h)) h$beta_ref else NA_real_
      er[[paste0("regenie_", rn, "_", s, "_se")]] <- if (nrow(h)) h$se else NA_real_
    }
    estimates[[length(estimates) + 1L]] <- er
  }
  rm(G); invisible(gc())
}
EX <- rbindlist(selected, fill = TRUE); ES <- rbindlist(estimates)
write_result(EX, "03_testosterone/recovery", "exemplar_selected")
write_result(rbindlist(profiles), "03_testosterone/recovery", "exemplar_profiles")
write_result(ES, "03_testosterone/recovery", "exemplar_estimates")

## 3. The GWAS of sex at every exemplar, aligned to the exemplar effect allele -------------
EX <- EX[!is.na(snp_id) & nzchar(snp_id), .(weight_design, category, snp_id, chr, bp)]
alleles <- unique(ES[arm == "testosterone", .(weight_design, snp_id, effect_allele, other_allele)])
EX <- merge(EX, alleles, by = c("weight_design", "snp_id"), all.x = TRUE, sort = FALSE)
EX[, allele_source := ifelse(is.na(effect_allele), NA_character_, "exemplar_estimates")]
# Exemplars chosen only under the other design take their alleles from its Stage 2 results.
for (i in which(is.na(EX$effect_allele))) {
  g <- fread(stage2_file(cfg, EX$weight_design[i], "testosterone", as.integer(EX$chr[i])),
             select = c("snp_id", "effect_allele", "other_allele"), showProgress = FALSE)[snp_id == EX$snp_id[i]]
  EX[i, `:=`(effect_allele = g$effect_allele, other_allele = g$other_allele,
             allele_source = paste0("stage2 ", EX$weight_design[i], "/testosterone"))]
}
S <- fread(p8$sex_reference, select = c("SNP", "A1", "A2", "N", "P", "BETA", "SE"), showProgress = FALSE)[SNP %in% EX$snp_id]
m <- match(EX$snp_id, S$SNP)
flip <- allele_sign(EX$effect_allele, EX$other_allele, S$A1[m], S$A2[m])
z <- S$BETA[m] / S$SE[m]
mlp <- ifelse(is.finite(S$P[m]) & S$P[m] > 0, -log10(S$P[m]), mlog10_p_chisq(z^2, 1))
OUT <- EX[, .(weight_design, category, snp_id, chr, bp, effect_allele, other_allele)]
OUT[, `:=`(sex_beta = S$BETA[m] * flip, sex_se = S$SE[m], sex_mlog10p = mlp, sex_n = as.numeric(S$N[m]),
           sex_z = z * flip, sex_allele_flip = flip, allele_source = EX$allele_source,
           gamma_sign_vs_sex = p8$gamma_sign_vs_sex)]
OUT[, sex_beta_gamma_oriented := sex_beta * gamma_sign_vs_sex]   # on the scale of gamma (male = -1)
OUT[is.na(flip), `:=`(sex_mlog10p = NA_real_, sex_n = NA_real_, sex_se = NA_real_)]
setorder(OUT, weight_design, category)
write_result(OUT, "03_testosterone/exemplars", "exemplar_sex_reference")
