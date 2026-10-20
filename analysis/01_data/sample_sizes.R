#!/usr/bin/env Rscript
# Sample sizes of every Stage 1 GWAS
#
# What: counts, for every trait with a Stage 1 quantile GWAS, the participants left at each
#   step of sample selection: keep list -> valid phenotype and complete covariates ->
#   present in the REGENIE leave-one-chromosome-out (LOCO) predictions (LOCO runs only) ->
#   matched to the genotype sample file = N analysed. Counts come from each chromosome's
#   Stage 1 log and covariance metadata; they must be identical on all 22 chromosomes.
# Paper: Supplementary Table S7
# Input: Stage 1 quantile GWAS logs and covariance metadata (one per trait and chromosome)
# Output: results/01_data/sample_sizes/{cohort_funnel,cohort_chr_counts}.csv
# Runs: on a cluster with the Stage 1 output, in a few minutes

.repo <- Sys.getenv("GMOQ_REPO", system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE))
source(file.path(.repo, "R", "load.R"))
cfg <- read_config()

## 1. Counts per trait and chromosome
traits <- stage1_traits(cfg)
per_chr <- rbindlist(lapply(traits$trait, function(tr) {
  files <- stage1_chr_paths(cfg, tr)
  rbindlist(lapply(seq_len(nrow(files)), function(j) {
    f <- files[j]
    log <- parse_stage1_log(f$stage1_log)
    meta <- if (file.exists(f$sigma0_meta)) jsonlite::read_json(f$sigma0_meta) else list()
    nuisance_n <- if (!is.na(f$rif_nuisance)) jsonlite::read_json(f$rif_nuisance)$n else NA
    data.table(
      trait = tr, chr = f$chr, source_kind = f$source_kind, log_found = log$log_found,
      n_keep = log$n_keep %||% NA, n_valid_pheno = log$n_valid_pheno %||% NA,
      n_absent_loco = log$n_absent_loco %||% NA, n_before_bgen = log$n_before_bgen %||% NA,
      n_bgen_matched = log$n_bgen_matched %||% NA,
      n_sigma0_meta = as.numeric(meta$n_samples %||% NA),
      n_rif_nuisance = as.numeric(nuisance_n %||% NA),
      n_progress = f$progress_n,
      n_covariates = as.numeric(meta$n_covariates %||% NA),
      covariates = paste(unlist(meta$covariates), collapse = " "),
      loco_mode = meta$loco_mode %||% "none",
      density_method = meta$density_method %||% NA_character_,
      pheno_col = meta$phenotype %||% NA_character_,
      pheno_file = log$pheno_file %||% NA, keep_file = log$keep_file %||% NA,
      covar_file = log$covar_file %||% NA, loco_file = log$loco_file %||% NA)
  }))
}), fill = TRUE)

## 2. One row per trait
# Every count should take one value across chromosomes; n_distinct_counts > 1 flags a
# trait where it does not.
single <- function(x) { v <- unique(x[!is.na(x)]); if (length(v) == 1L) v else x[NA_integer_][1] }
n_distinct <- function(x) uniqueN(x[!is.na(x)])
funnel <- per_chr[, .(
  n_chr = .N, n_logs = sum(log_found),
  n_keep = single(n_keep), n_valid_pheno = single(n_valid_pheno),
  n_absent_loco = single(n_absent_loco), n_before_bgen = single(n_before_bgen),
  n_analysed = single(n_bgen_matched), n_sigma0_meta = single(n_sigma0_meta),
  n_rif_nuisance = single(n_rif_nuisance),
  n_covariates = single(n_covariates), covariates = single(covariates),
  loco_mode = single(loco_mode), density_method = single(density_method),
  pheno_col = single(pheno_col), pheno_file = single(pheno_file), keep_file = single(keep_file),
  covar_file = single(covar_file), loco_file_chr1 = loco_file[chr == min(chr)][1],
  n_distinct_counts = max(n_distinct(n_keep), n_distinct(n_valid_pheno),
                          n_distinct(n_bgen_matched), n_distinct(n_sigma0_meta)),
  n_distinct_covariates = n_distinct(covariates),
  source_kinds = paste(sort(unique(source_kind)), collapse = ";")), by = trait]
funnel <- merge(traits[, .(trait, group, n_config, fig6)], funnel, by = "trait", sort = FALSE)

# N_phenotyped (Table S7) is n_valid_pheno; N analysed is the genotype-matched count, which
# must equal the N recorded with the pooled covariance.
funnel[, n_after_loco := ifelse(is.na(n_absent_loco), n_valid_pheno, n_valid_pheno - n_absent_loco)]
funnel[, n_dropped_bgen := n_before_bgen - n_analysed]
funnel[, analysed_equals_sigma0 := n_analysed == n_sigma0_meta]
funnel[, analysed_equals_config := n_analysed == n_config]
funnel[, consistent_across_chr := n_distinct_counts == 1L & n_distinct_covariates == 1L & n_logs == n_chr]
setcolorder(funnel, c("trait", "group", "fig6", "n_keep", "n_valid_pheno", "n_absent_loco",
                      "n_after_loco", "n_before_bgen", "n_dropped_bgen", "n_analysed",
                      "n_sigma0_meta", "n_rif_nuisance", "n_config"))
print(funnel[, .(trait, n_keep, n_valid_pheno, n_after_loco, n_analysed, n_config,
                 n_covariates, loco_mode, consistent_across_chr)])

write_result(funnel, "01_data/sample_sizes", "cohort_funnel")
write_result(per_chr[, !c("pheno_file", "keep_file", "covar_file", "loco_file")],
             "01_data/sample_sizes", "cohort_chr_counts")
