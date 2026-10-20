#!/usr/bin/env Rscript
# Fit the two-component skew-normal mixtures and build their Stage 2 weights
#
# What: fits a two-component skew-normal mixture to each of six phenotypes (log testosterone,
#   unshrunk and with the male values shrunk towards the female mean by 0.2, 0.5 and 0.7; log
#   JCV VP1 and log MCV VP1 antibody levels) on exactly the participants Stage 1 analysed. From
#   each fit it builds the Stage 2 weight matrix for the three mixture parameters (membership
#   gamma, lower-component location beta_left, upper-component location beta_right) under the
#   two weight designs compared in Supplementary Note 2: 'model' divides every column by the
#   fitted mixture density; 'hybrid' divides the membership column by the Stage 1 kernel density.
#   gamma is logit P(lower component); the paper reports upper-component membership, -gamma.
# Paper: Fig. 1b; Fig. 2; Supplementary Fig. S7a; Supplementary Note 2; Methods
# Input: UK Biobank phenotype, covariate and sample files of each Stage 1 GWAS; Stage 1 output
#   per chromosome (rif_nuisance.json: quantile levels, baseline quantiles and kernel density)
# Output: results/02_mixtures/stage2_mixtures/{weight_matrices,weight_density}.csv;
#   per-chromosome weight matrices and the fitted mixtures in the working directory kept with the UK Biobank data
# Runs: on a cluster with UK Biobank data (about 1 hour; the mixture fit to 399,710 values
#   dominates)

repo <- Sys.getenv("GMOQ_REPO", system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE))
source(file.path(repo, "R", "load.R"))
cfg <- read_config()
suppressPackageStartupMessages(library(gmoq))   # the paper used a development version of this package (Stage 2 commit 9fee67d)
taus_grid <- round(tau_grid(cfg), 3)

weight_rows <- list(); density_rows <- list()
for (i in seq_len(nrow(MIXTURE_FITS))) {
  f <- MIXTURE_FITS[i]

  ## 1. Fit the mixture to the analysed sample ----------------------------------------
  # The sample, and therefore the fit, is the same for every chromosome (LOCO predictions
  # enter Stage 1 as a covariate, not by changing the phenotype), so the mixture is fitted
  # once per phenotype.
  inp <- stage1_inputs(cfg, f$fit)
  sample <- stage1_analysed_sample(inp$pheno, inp$trait, inp$covar, inp$keep, inp$loco, inp$bgen_sample)
  y <- sample$y
  # The rebuilt sample must be the one Stage 1 used: same size, mean and standard deviation.
  s1 <- read_json(file.path(stage1_dir(cfg, f$cohort, 1), "rif_nuisance.json"))
  stopifnot(length(y) == s1$n, isTRUE(all.equal(mean(y), s1$mean, tolerance = 1e-8)),
            isTRUE(all.equal(sd(y), s1$sd, tolerance = 1e-8)))
  components <- fit_skew_normal_mixture(y, seed = MIXTURE_SEED)
  saveRDS(list(fit = f$fit, components = components, n = length(y), seed = MIXTURE_SEED),
          mixture_fit_file(cfg, f$fit))
  message(f$fit, ": mixture fitted to ", length(y), " participants")

  ## 2. Weight matrices under both designs, per chromosome ----------------------------
  # Stage 2 uses the quantile levels, baseline quantiles q_tau and kernel density that Stage 1
  # used on each chromosome (identical across chromosomes up to rounding).
  for (design in WEIGHT_DESIGNS) {
    W_by_chr <- list()
    for (chr in 1:22) {
      dir <- stage1_dir(cfg, f$cohort, chr)
      density_file <- file.path(dir, "rif_nuisance.json")
      rn <- read_json(density_file)
      taus <- as.numeric(rn$taus); q_tau <- as.numeric(rn$q_tau); f_kde <- as.numeric(rn$fhat)
      stopifnot(length(taus) == 45L, max(abs(round(taus, 3) - taus_grid)) <= 5e-7)
      W <- mixture_weight_matrix(components, taus, q_tau, design, density_file)
      if (!check_rank(unclass(W))$full) stop("W is not of full column rank: ", f$fit, " chr", chr)
      saveRDS(list(W = W, taus = round(taus, 3), stage1_dir = dir, components = components),
              weight_file(cfg, design, f$fit, chr))
      W_by_chr[[chr]] <- W
      # Fitted mixture density against the Stage 1 kernel density at each q_tau (Fig. S7a);
      # it does not depend on the design.
      mf <- mixture_functions(components)
      if (design == "model") density_rows[[length(density_rows) + 1L]] <- data.table(
        fit = f$fit, chr = chr, tau = round(taus, 3), q_tau = q_tau, f_kde = f_kde,
        f_model = pmax(mf$pdf(q_tau, mf$p), 1e-12), F_model_minus_tau = mf$cdf(q_tau, mf$p) - taus)
    }

    ## 3. Weight matrix summary across chromosomes -------------------------------------
    W1 <- unclass(W_by_chr[[1]])
    W_all <- simplify2array(lapply(W_by_chr, function(w) matrix(unclass(w), nrow(w), ncol(w))))
    for (j in seq_len(ncol(W1))) weight_rows[[length(weight_rows) + 1L]] <- data.table(
      fit = f$fit, weight_design = design, tau = taus_grid, param = colnames(W1)[j],
      W_chr01 = as.numeric(W1[, j]),
      W_min_chr = apply(W_all[, j, , drop = FALSE], 1, min),
      W_max_chr = apply(W_all[, j, , drop = FALSE], 1, max))
  }
}

## 4. Write the aggregate tables ------------------------------------------------------
write_result(rbindlist(weight_rows), "02_mixtures/stage2_mixtures", "weight_matrices")

# Density table: chromosome 1 values and the range of the density ratio over chromosomes.
D <- rbindlist(density_rows)
D[, ratio_model_over_kde := f_model / f_kde]
over_chr <- D[, .(ratio_min_chr = min(ratio_model_over_kde), ratio_max_chr = max(ratio_model_over_kde),
                  f_kde_rel_range_chr = (max(f_kde) - min(f_kde)) / mean(f_kde)), by = .(fit, tau)]
write_result(merge(D[chr == 1, !"chr"], over_chr, by = c("fit", "tau")),
             "02_mixtures/stage2_mixtures", "weight_density")
