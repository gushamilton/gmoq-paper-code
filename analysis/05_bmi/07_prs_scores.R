#!/usr/bin/env Rscript
# Leave-one-chromosome-out BMI polygenic scores for the interaction analysis
#
# What: builds the participant-level sample for the BMI x polygenic score analysis: the same
#   461,320 participants as the BMI Stage 1 GWAS, with BMI, log BMI, age, sex and ten genetic
#   principal components from the Stage 1 input files. The polygenic score for a lead on
#   chromosome c is the sum of the clumping-and-thresholding partial scores (analysis/00_gwas/06_bmi_prs_clump_threshold.sh)
#   of the other 21 chromosomes,
#   standardised over the sample, so the tested variant is never in LD with its own score.
#   The sample is also split at random into halves for the cross-fitted link test.
# Paper: Fig. 6; Supplementary Fig. S14; Methods (BMI x PRS interaction)
# Input: Stage 1 BMI sample (keep list), phenotype and covariate files; partial scores per
#   chromosome; lead variant list (results/01_data/lead_profiles)
# Output: participant-level files in the private area only (analysis sample with the 22
#   scores and the half label; the same without identifiers for the RIF step; keep lists of
#   the halves; lead rsids per chromosome). Nothing is written to results/.
# Runs: on a cluster with UK Biobank data, in a few minutes

repo <- Sys.getenv("GMOQ_REPO", unset = NA)
if (is.na(repo)) repo <- system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE)
source(file.path(repo, "R", "load.R"))
cfg <- read_config()
P <- cfg$prs
PCS <- paste0("PC", 1:10)
private <- file.path(cfg$bp$private_root, "bmi_prs", "prep")
dir.create(file.path(private, "snps"), recursive = TRUE, showWarnings = FALSE)

## 1. Analysis sample ------------------------------------------------------------------------
keep <- fread(P$keep_matched, header = FALSE, col.names = c("FID", "IID"), colClasses = "character")
pheno <- fread(P$pheno, select = c("FID", "IID", "bmi", "log_bmi"), colClasses = list(character = c("FID", "IID")))
covar <- fread(P$covar, colClasses = list(character = c("FID", "IID")))
A <- merge(keep, pheno, by = c("FID", "IID"), sort = FALSE)
A <- merge(A, covar[, c("FID", "IID", "age", "sex", PCS), with = FALSE], by = c("FID", "IID"), sort = FALSE)
A <- A[match(keep$IID, A$IID)]                               # keep-file order
n <- nrow(A)
stopifnot(n == P$n, !anyNA(A))

## 2. Leave-one-chromosome-out scores ---------------------------------------------------------
partial <- matrix(NA_real_, n, 22)
for (c in 1:22) {
  sc <- fread(file.path(P$partial_scores, sprintf("chr%02d_partial.sscore", c)), select = c("IID", "BETA_SUM"),
              colClasses = list(character = "IID"))
  i <- match(A$IID, sc$IID)
  partial[, c] <- ifelse(is.na(i), 0, sc$BETA_SUM[i])
}
total <- rowSums(partial)
for (c in 1:22) {
  s <- total - partial[, c]
  set(A, j = sprintf("PRS%02d", c), value = (s - mean(s)) / stats::sd(s))
}

## 3. Random halves (cross-fitted link test) ---------------------------------------------------
set.seed(P$link_test_seed)
A[, half := sample(rep(1:2, length.out = n))]

## 4. Private outputs ----------------------------------------------------------------------
saveRDS(A, file.path(private, "analysis.rds"))
fwrite(A[, c("bmi", "log_bmi", "half", "age", "sex", PCS, sprintf("PRS%02d", 1:22)), with = FALSE],
       file.path(private, "rif_input.csv.gz"))
for (h in 1:2) fwrite(A[half == h, .(FID, IID)], file.path(private, sprintf("keep_half%d.txt", h)), sep = " ", col.names = FALSE)
leads <- as.data.table(read_result("01_data/lead_profiles", "lead_profiles"))[lead_set == P$lead_set]
for (c in 1:22) writeLines(leads[chr == c]$snp_id, file.path(private, "snps", sprintf("chr%02d.txt", c)))
