#!/usr/bin/env Rscript
# Distribution of log testosterone, overall and by sex
#
# What: histograms of the analysed testosterone phenotype (log scale) for the 399,710
#   participants of the Stage 1 GWAS, unshrunk and in the three arms in which male values were
#   shrunk towards the female mean by F = 0.2, 0.5 and 0.7, overall and separately in women and
#   men. These show how far the two sexes, the two mixture components, overlap in each arm.
#   Bins are merged until every bin holds at least 10 participants.
# Paper: Fig. 1a
# Input: UK Biobank testosterone phenotype file, covariates, sex, keep list, REGENIE LOCO
#   predictions, imputed-genotype sample file (paths in config.yaml, section p08)
# Output: results/03_testosterone/recovery/phenotype_histograms.csv
# Runs: on a cluster with UK Biobank data, in about 2 minutes

repo <- Sys.getenv("GMOQ_REPO", system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE))
source(file.path(repo, "R", "load.R"))
cfg <- read_config()
p8 <- cfg$p08
arms <- c("testosterone", "testosterone_shrink_0.2", "testosterone_shrink_0.5", "testosterone_shrink_0.7")
sex <- fread(p8$sex_covar_file, select = c("FID", "sex"), colClasses = list(character = "FID"))

rows <- list()
for (arm in arms) {
  ## 1. The analysed sample (the Stage 1 sample rule) ------------------------------------
  s <- stage1_analysed_sample(p8$phenotype_file, arm, p8$covar_file, p8$keep_file, p8$loco_list, p8$bgen_sample)
  if (nrow(s) != cfg$cohorts$testosterone$n) stop(arm, ": analysed sample has ", nrow(s), " participants")
  s[, sex := sex$sex[match(FID, sex$FID)]]           # 0 = female, 1 = male

  ## 2. Histograms: inner edges between the 0.1st and 99.9th percentiles ------------------
  edges <- as.numeric(stats::quantile(s$y, c(0.001, 0.999)))
  for (grp in c("all", "female", "male")) {
    y <- switch(grp, all = s$y, female = s$y[s$sex == 0], male = s$y[s$sex == 1])
    h <- privacy_histogram(y, edges[1], edges[2], p8$histogram_bins, cfg$privacy$min_cell_count)
    rows[[length(rows) + 1L]] <- data.table(phenotype = arm, group = grp, h[, .(bin_lower, bin_upper, count)],
                                            n_group = length(y), density = h$density)
  }
}
write_result(rbindlist(rows), "03_testosterone/recovery", "phenotype_histograms")
