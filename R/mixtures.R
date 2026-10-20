# Shared functions for the two-component mixture analyses
#
# What: functions used by the mixture fits (02_mixtures), the testosterone analyses
#   (03_testosterone) and the antibody analyses (04_antibodies): which traits are fitted, where
#   their Stage 1 output lives, how the analysed sample is rebuilt, how the skew-normal mixture
#   is fitted, and how its Stage 2 weight matrix is built under the two weight designs.
# Paper: Methods (two-component mixtures); Supplementary Note 2
# Loaded by R/load.R.
#
# Sign convention used throughout the code. The mixture parameters are
#   gamma       = logit P(lower component)  (per-allele effect on the log-odds of membership
#                 of the LOWER component: the female component for testosterone, the
#                 seronegative-like component for the antibodies),
#   beta_left   = shift of the lower component's location,
#   beta_right  = shift of the upper component's location.
# The paper reports membership of the UPPER component, which is -gamma. Output files keep the
# code's gamma; figures and tables flip the sign where they show membership.

suppressPackageStartupMessages(library(data.table))

## 1. The six mixture fits ----------------------------------------------------------------
# fit:    name used in file names and results;
# cohort: key of config.yaml cohorts (location of the Stage 1 GWAS);
# trait:  phenotype column in the phenotype file.
MIXTURE_FITS <- data.table(
  fit    = c("testosterone", "testosterone_shrink_0_2", "testosterone_shrink_0_5",
             "testosterone_shrink_0_7", "jcv_vp1_log_all", "mcv_vp1_log_all"),
  cohort = c("testosterone", "testosterone_shrink_0_2", "testosterone_shrink_0_5",
             "testosterone_shrink_0_7", "jcv_vp1_log_all", "mcv_vp1_log_all"),
  trait  = c("testosterone", "testosterone_shrink_0.2", "testosterone_shrink_0.5",
             "testosterone_shrink_0.7", "jcv_vp1_log_all", "mcv_vp1_log_all"),
  model  = c(rep("skew", 4), "mixture", "mixture"))
MIXTURE_PARAMS <- c("gamma", "beta_left", "beta_right")
WEIGHT_DESIGNS <- c("hybrid", "model")
MIXTURE_SEED <- 42L      # seed of the mixture EM initialisation (k-means start values)

## 2. Locations ---------------------------------------------------------------------------
# Directory of the Stage 1 output for one cohort and chromosome. Unshrunk testosterone
# chromosomes 1 and 22 were produced in an earlier run with the same settings and sample, so
# they are read from a second location when absent from the first.
stage1_dir <- function(cfg, cohort, chr) {
  d <- file.path(cfg$cohorts[[cohort]]$root, sprintf("chr%02d", chr))
  if (file.exists(file.path(d, "scan.stage1.tsv.gz"))) return(d)
  d2 <- file.path(cfg$p07$stage1_fallback_root, cohort, sprintf("chr%02d", chr))
  if (file.exists(file.path(d2, "scan.stage1.tsv.gz"))) return(d2)
  stop("no Stage 1 output for ", cohort, " chromosome ", chr)
}

# Working directory for per-variant and participant-level intermediates. These stay with the
# UK Biobank data and are never written to results/.
private_dir <- function(cfg, ...) {
  d <- file.path(cfg$bp$private_root, "gmoq-paper-code", ...)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}
mixture_fit_file <- function(cfg, fit) file.path(private_dir(cfg, "02_mixtures", "mixture_fits"),
                                                 paste0(fit, ".rds"))
weight_file <- function(cfg, design, fit, chr)
  file.path(private_dir(cfg, "02_mixtures", "weights", design, fit), sprintf("chr%02d.rds", chr))
stage2_file <- function(cfg, design, fit, chr)
  file.path(private_dir(cfg, "02_mixtures", "stage2", design, fit), sprintf("chr%02d.tsv.gz", chr))

# Per-variant Stage 2 results of one fit and design, all chromosomes.
read_stage2 <- function(cfg, design, fit, select = NULL) {
  rbindlist(lapply(1:22, function(chr)
    fread(stage2_file(cfg, design, fit, chr), select = select, showProgress = FALSE)))
}

read_json <- function(path) jsonlite::read_json(path, simplifyVector = TRUE)
read_stage1 <- function(dir) fread(cmd = paste("gzip -dc", shQuote(file.path(dir, "scan.stage1.tsv.gz"))),
                                   showProgress = FALSE)

## 3. The analysed sample -----------------------------------------------------------------
# Rebuilds the sample that Stage 1 analysed, in phenotype-file order: participants in the
# keep list with complete covariates and a finite phenotype, present in the REGENIE
# leave-one-chromosome-out (LOCO) predictions and in the imputed-genotype sample file.
#   loco: a REGENIE .loco file or a *_pred.list naming them (membership only is used; every
#         .loco file of one step 1 run has the same participants).
# Value: data.table FID, IID, y.
stage1_analysed_sample <- function(pheno_file, trait, covar_file, keep_file, loco, bgen_sample) {
  keep <- unique(fread(keep_file, header = FALSE, colClasses = "character")[[1]])
  ph <- fread(pheno_file, select = c("FID", "IID", trait), colClasses = list(character = c("FID", "IID")))
  setnames(ph, trait, "y")
  cv <- fread(covar_file, colClasses = list(character = c("FID", "IID")))
  ph <- ph[FID %in% keep & FID %in% cv$FID]
  cvm <- cv[match(ph$FID, cv$FID)]
  complete_cov <- Reduce(`&`, lapply(setdiff(names(cv), c("FID", "IID")),
                                     function(v) is.finite(as.numeric(cvm[[v]]))))
  ph <- ph[complete_cov & is.finite(as.numeric(y))]
  loco_file <- if (grepl("\\.loco$", loco)) loco else strsplit(readLines(loco)[1], "\\s+")[[1]][2]
  loco_ids <- scan(text = readLines(loco_file, n = 1L), what = character(), quiet = TRUE)[-1]
  # REGENIE may label participants by FID, IID or a combination of the two.
  in_loco <- ph$FID %in% loco_ids | ph$IID %in% loco_ids |
    paste(ph$FID, ph$IID, sep = "_") %in% loco_ids | paste(ph$FID, ph$IID, sep = " ") %in% loco_ids |
    paste(ph$FID, ph$IID, sep = ":") %in% loco_ids
  bgen_ids <- fread(bgen_sample, skip = 2L, header = FALSE, colClasses = "character")[[1]]
  ph <- ph[in_loco & FID %in% bgen_ids]
  ph[, y := as.numeric(y)]
  ph[]
}

# Input files of the Stage 1 GWAS of one fit. Testosterone uses the paper phenotype files and
# the REGENIE step 1 of its own arm; the antibody Stage 1 recorded its inputs in inputs.log.
stage1_inputs <- function(cfg, fit) {
  f <- MIXTURE_FITS[match(fit, MIXTURE_FITS$fit)]
  if (startsWith(f$fit, "testosterone")) {
    step1_root <- dirname(dirname(cfg$p08$loco_list))
    return(list(pheno = cfg$p08$phenotype_file, covar = cfg$p08$covar_file, keep = cfg$p08$keep_file,
                loco = file.path(step1_root, f$fit, paste0(f$fit, "_step1_pred.list")),
                bgen_sample = cfg$p08$bgen_sample, trait = f$trait))
  }
  log <- readLines(file.path(stage1_dir(cfg, f$cohort, 22), "inputs.log"))
  from <- function(pattern) sub(".* from ([^ ]+).*$", "\\1", grep(pattern, log, value = TRUE)[1])
  list(pheno = from("Loading phenotypes from"), covar = from("Loading covariates from"),
       keep = from("Loading sample keep list from"), loco = from("Loading LOCO predictions for"),
       bgen_sample = cfg$p09$bgen_sample, trait = f$trait)
}

## 4. The skew-normal mixture -------------------------------------------------------------
# Two-component skew-normal mixture fitted by EM (mixsmsn::smsn.mix) to the phenotype of the
# analysed sample, with components ordered by location so that "lower" and "upper" mean the
# same thing on every run.
# Value: data.table component (lower, upper), weight, xi (location), omega (scale), alpha (shape).
fit_skew_normal_mixture <- function(y, seed = MIXTURE_SEED) {
  set.seed(seed)
  fit <- mixsmsn::smsn.mix(y = y, g = 2, family = "Skew.normal", nu = 1, get.init = TRUE,
                           iter.max = 500, error = 1e-4, calc.im = FALSE)
  ord <- order(as.numeric(fit$mu))
  comp <- data.table(component = c("lower", "upper"), weight = as.numeric(fit$pii)[ord],
                     xi = as.numeric(fit$mu)[ord], omega = sqrt(as.numeric(fit$sigma2)[ord]),
                     alpha = as.numeric(fit$shape)[ord])
  if (any(!is.finite(unlist(comp[, -1]))) || any(comp$weight <= 0))
    stop("the skew-normal mixture fit returned invalid parameters")
  comp
}

# Distribution function, density and the derivatives of the distribution function with respect
# to the three Stage 2 parameters, for make_weights_generic(). p1 is P(lower component);
# because gamma = logit(p1), dF/dgamma = dF/dp1 * p1 (1 - p1). The location derivatives are
# central finite differences of F in xi1 and xi2.
mixture_functions <- function(components) {
  cm <- as.data.frame(components)
  lo <- cm[cm$component == "lower", ]; up <- cm[cm$component == "upper", ]
  p <- list(p1 = lo$weight, xi1 = lo$xi, omega1 = lo$omega, alpha1 = lo$alpha,
            xi2 = up$xi, omega2 = up$omega, alpha2 = up$alpha)
  cdf <- function(x, z) z$p1 * sn::psn(x, xi = z$xi1, omega = z$omega1, alpha = z$alpha1) +
    (1 - z$p1) * sn::psn(x, xi = z$xi2, omega = z$omega2, alpha = z$alpha2)
  pdf <- function(x, z) z$p1 * sn::dsn(x, xi = z$xi1, omega = z$omega1, alpha = z$alpha1) +
    (1 - z$p1) * sn::dsn(x, xi = z$xi2, omega = z$omega2, alpha = z$alpha2)
  dF_dp1 <- gmoq::make_fd_grad("p1", cdf)
  list(p = p, cdf = cdf, pdf = pdf, grads = list(
    gamma = function(x, z) dF_dp1(x, z) * z$p1 * (1 - z$p1),
    beta_left = gmoq::make_fd_grad("xi1", cdf),
    beta_right = gmoq::make_fd_grad("xi2", cdf)))
}

# Stage 2 weight matrix W (quantile levels x {gamma, beta_left, beta_right}):
#   W_j(tau) = -dF(q_tau)/dtheta_j / f(q_tau).
# design = "model":  every column is divided by the fitted mixture density f_model.
# design = "hybrid": the membership column is divided by the Stage 1 kernel density instead
#   (stage1_density_file = Stage 1 rif_nuisance.json, which holds the density Stage 1 divided
#   by when it built the recentred influence function). Stage 1 estimates each quantile effect
#   as a distribution-function slope divided by that kernel density, so the membership column,
#   which is a difference of component distribution functions, has the scale of the Stage 1
#   output only if it is divided by the same density. The location columns keep f_model: they
#   are then the component shares P(k | q_tau), which sum to 1 at every tau.
mixture_weight_matrix <- function(components, taus, q_tau, design, stage1_density_file = NULL) {
  mf <- mixture_functions(components)
  gmoq::make_weights_generic(
    taus = taus, q_tau = q_tau, dist_cdf = mf$cdf, dist_pdf = mf$pdf, params = mf$p,
    grad_funcs = mf$grads, density = design,
    stage1_density = if (design == "model") NULL else stage1_density_file,
    stage1_columns = "gamma", unit_sum = list(c("beta_left", "beta_right")), identity_tol = 1e-8)
}

# Skew-normal density (as sn::dsn), used for plotting the fitted components.
dsn_density <- function(x, xi, omega, alpha) 2 / omega * dnorm((x - xi) / omega) * pnorm(alpha * (x - xi) / omega)

## 5. Allele alignment ---------------------------------------------------------------------
# +1 when (effect, other) equals (a1, a0), -1 when swapped, NA otherwise (no strand flips).
allele_sign <- function(effect, other, a1, a0) {
  s <- rep(NA_real_, length(effect))
  s[effect == a1 & other == a0] <- 1
  s[effect == a0 & other == a1] <- -1
  s
}

## 6. Comparison of a G-MoQ parameter with a reference GWAS -------------------------------
# Joins a G-MoQ parameter (multiplied by sign) to a reference GWAS (snp_id, a1 = effect allele,
# a2, n, beta, se, mlog10p) on variant identifier and aligns the reference to the G-MoQ effect
# allele; variants whose alleles do not match either way are dropped.
# Value: one row per variant with b, s (G-MoQ), beta_ref, se (reference), z_g, z_r.
harmonise_with_reference <- function(g, param, sign, ref) {
  m <- merge(g[, .(snp_id, chr, bp, ea = effect_allele, oa = other_allele,
                   b = get(param) * sign, s = get(paste0(param, "_se")))], ref, by = "snp_id")
  flip <- ifelse(m$ea == m$a1 & m$oa == m$a2, 1, ifelse(m$ea == m$a2 & m$oa == m$a1, -1, NA))
  m[, beta_ref := beta * flip]
  m <- m[!is.na(flip) & is.finite(b) & is.finite(beta_ref) & se > 0 & s > 0]
  m[, `:=`(z_g = b / s, z_r = beta_ref / se)]
  m
}
lambda_from_z <- function(z) stats::median(z^2, na.rm = TRUE) / qchisq(0.5, 1)

# Agreement of the harmonised effects: Pearson correlation of effects and of z statistics,
# least-squares slope of G-MoQ on the reference, median ratio of standard errors, and the
# genomic inflation of both on the common variants.
concordance_stats <- function(m) {
  if (nrow(m) < 3) return(data.table(n = nrow(m)))
  f <- lm(b ~ beta_ref, m)
  data.table(n = nrow(m), pearson_beta = cor(m$b, m$beta_ref), pearson_z = cor(m$z_g, m$z_r),
             slope_gmoq_on_regenie = coef(f)[[2]], slope_se = summary(f)$coefficients[2, 2],
             intercept = coef(f)[[1]], median_se_ratio_gmoq_over_regenie = median(m$s / m$se),
             lambda_gc_gmoq_common = lambda_from_z(m$z_g), lambda_gc_regenie_common = lambda_from_z(m$z_r))
}

## 7. LD score regression ------------------------------------------------------------------
# munge_sumstats.py on a summary-statistics file (SNP, A1, A2, N, P, BETA), HapMap3 alleles.
ldsc_munge <- function(cfg, raw_file, out_prefix) run_shell(sprintf(
  "module load %s && munge_sumstats.py --sumstats %s --out %s --merge-alleles %s --snp SNP --a1 A1 --a2 A2 --p P --N-col N --signed-sumstats BETA,0 --ignore SE",
  cfg$p08$ldsc_module, raw_file, out_prefix, cfg$p08$ldsc_merge_alleles))

# Genetic correlation of two munged files (1000 Genomes European LD scores); the log is
# <out_prefix>.log. An existing complete log is reused.
ldsc_rg <- function(cfg, munged_1, munged_2, out_prefix) {
  log <- paste0(out_prefix, ".log")
  if (!file.exists(log) || !length(grep("^\\s*p1\\s+p2\\s+rg", readLines(log))))
    run_shell(sprintf("module load %s && ldsc.py --rg %s,%s --ref-ld-chr %s --w-ld-chr %s --out %s",
                      cfg$p08$ldsc_module, munged_1, munged_2, cfg$p08$ldsc_ld_dir, cfg$p08$ldsc_ld_dir, out_prefix))
  readLines(log)
}

# The rg summary line of an ldsc.py --rg log, as a named list (rg, se, z, p, h2_obs, ...).
ldsc_rg_summary <- function(log) {
  i <- grep("^\\s*p1\\s+p2\\s+rg", log)
  if (!length(i)) stop("no rg summary in the LDSC log")
  setNames(as.list(strsplit(trimws(log[i[1] + 1]), "\\s+")[[1]]), strsplit(trimws(log[i[1]]), "\\s+")[[1]])
}
ldsc_n_snps <- function(log) suppressWarnings(as.integer(sub(".*?([0-9]+) SNPs with valid alleles.*", "\\1",
                                                         grep("SNPs with valid alleles", log, value = TRUE)[1])))

# Heritability block of phenotype 1 or 2 from an ldsc.py --rg log: observed-scale h2, lambda
# GC, mean chi2, intercept and ratio (intercept - 1) / (mean chi2 - 1).
ldsc_h2_block <- function(log, which) {
  i0 <- grep(if (which == 1) "^Heritability of phenotype 1\\s*$" else "^Heritability of phenotype 2", log)[1]
  if (is.na(i0)) return(NULL)
  blk <- log[i0:min(length(log), i0 + 8)]
  num <- function(pat) {
    x <- grep(pat, blk, value = TRUE)[1]
    if (is.na(x)) return(c(NA_real_, NA_real_))
    x <- sub("^[^:]*:", "", x)
    v <- regmatches(x, gregexpr("-?[0-9.]+(e-?[0-9]+)?", x))[[1]]
    c(as.numeric(v[1]), if (length(v) > 1) as.numeric(v[2]) else NA_real_)
  }
  h2 <- num("^Total Observed scale h2"); lg <- num("^Lambda GC"); mc <- num("^Mean Chi")
  ic <- num("^Intercept"); rt <- num("^Ratio:")
  data.table(h2 = h2[1], h2_se = h2[2], ldsc_lambda_gc = lg[1], mean_chi2 = mc[1],
             intercept = ic[1], intercept_se = ic[2], ratio = rt[1], ratio_se = rt[2])
}

## 8. Phenotype histograms that disclose no small counts ---------------------------------
# Inner bin edges are evenly spaced between two central quantiles; the two outer bins are
# open-ended, so no individual minimum or maximum is shown. Any bin with fewer than min_count
# participants is merged with its smaller neighbour until none remains.
# Value: data.table bin_lower, bin_upper, count, density (NA for the open outer bins).
privacy_histogram <- function(y, inner_from, inner_to, n_bins, min_count = 10L) {
  breaks <- c(-Inf, seq(inner_from, inner_to, length.out = n_bins + 1L), Inf)
  repeat {
    count <- tabulate(findInterval(y, breaks, rightmost.closed = TRUE, all.inside = TRUE), length(breaks) - 1L)
    small <- which(count < min_count)
    if (!length(small) || length(breaks) <= 3L) break
    j <- small[1]
    drop <- if (j == 1L) 2L else if (j == length(count)) length(breaks) - 1L else
      if (count[j - 1L] <= count[j + 1L]) j else j + 1L
    breaks <- breaks[-min(max(drop, 2L), length(breaks) - 1L)]
  }
  data.table(bin_lower = signif(head(breaks, -1), 7), bin_upper = signif(tail(breaks, -1), 7), count = count,
             density = ifelse(is.finite(diff(breaks)), count / (length(y) * diff(breaks)), NA_real_))
}

# Run a shell command through a login shell (so that cluster modules can be loaded); stop
# with the last lines of output if it fails.
run_shell <- function(cmd) {
  out <- system2("bash", c("-lc", shQuote(cmd)), stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status")
  if (!is.null(status) && status != 0)
    stop("command failed: ", cmd, "\n", paste(tail(out, 20), collapse = "\n"))
  out
}
