# Summaries of genome-wide results, and readers for the exported lead profiles.
#
# Per-variant genome-wide results stay with the data. What is exported is aggregate:
# distance-clumped lead rows, Manhattan points thinned to the significant variants plus a
# maximum per 1 Mb bin, QQ-plot quantiles and the genomic-control lambda.

# -log10 p of a chi-square statistic, computed on the log scale so that p < 1e-308 does not
# underflow.
mlog10_p_chisq <- function(stat, df) {
  -stats::pchisq(stat, df, lower.tail = FALSE, log.p = TRUE) / log(10)
}

# -log10 p of a Wald statistic with the Hotelling reference
# T (nu - k + 1) / (k nu) ~ F(k, nu - k + 1) (see hotelling_p()), on the log scale.
mlog10_p_hotelling <- function(T, k, nu) {
  -stats::pf(T * (nu - k + 1) / (k * nu), k, nu - k + 1, lower.tail = FALSE,
             log.p = TRUE) / log(10)
}

# Genomic-control lambda: the median chi-square implied by the p-values divided by the
# median of chi2(df). p-values outside (0, 1] are ignored.
lambda_gc <- function(p, df = 1) {
  p <- p[is.finite(p) & p > 0 & p <= 1]
  if (!length(p)) return(NA_real_)
  stats::qchisq(stats::median(p), df, lower.tail = FALSE) / stats::qchisq(0.5, df)
}

# Points for a QQ plot from -log10 p-values: the max_tail most significant ranks, plus
# n_points ranks spaced evenly on the log scale over the rest. Expected -log10 p at rank i
# of n is -log10((i - 0.5) / n).
# Value: data.frame(rank, expected_mlog10p, observed_mlog10p, n_snps).
qq_quantiles <- function(mlog10p, n_points = 400L, max_tail = 200L) {
  x <- sort(mlog10p[is.finite(mlog10p)], decreasing = TRUE)
  n <- length(x)
  if (!n) return(data.frame(rank = integer(0), expected_mlog10p = numeric(0),
                            observed_mlog10p = numeric(0), n_snps = integer(0)))
  idx <- unique(c(seq_len(min(max_tail, n)),
                  round(exp(seq(log(min(max_tail, n)), log(n), length.out = n_points)))))
  idx <- sort(unique(pmin(pmax(idx, 1L), n)))
  data.frame(rank = idx, expected_mlog10p = -log10((idx - 0.5) / n),
             observed_mlog10p = x[idx], n_snps = n)
}

# Thinned Manhattan plot data.
#   dt: data.frame with snp_id, chr, bp and the -log10 p columns named in cols;
#   keep_p: variants with p below this in any column are kept individually;
#   bin_bp: bin width for the per-chromosome bin maxima of every column.
# Value: list(points = kept variant rows, bins = chr, bin_start_bp, n_variants and
#   max_<col> per bin, keep_p, bin_bp).
manhattan_thin <- function(dt, cols, keep_p = 1e-3, bin_bp = 1e6) {
  dt <- as.data.frame(dt)
  thr <- -log10(keep_p)
  mx <- do.call(pmax, c(lapply(cols, function(cc) dt[[cc]]), na.rm = TRUE))
  keep <- is.finite(mx) & mx >= thr
  points <- dt[keep, c("snp_id", "chr", "bp", cols), drop = FALSE]
  bin <- floor(dt$bp / bin_bp)
  key <- paste(dt$chr, bin)
  agg <- lapply(cols, function(cc) tapply(dt[[cc]], key, max, na.rm = TRUE))
  k <- names(agg[[1]])
  sp <- do.call(rbind, strsplit(k, " ", fixed = TRUE))
  bins <- data.frame(chr = as.integer(sp[, 1]), bin_start_bp = as.numeric(sp[, 2]) * bin_bp,
                     n_variants = as.integer(table(key)[k]))
  for (i in seq_along(cols)) bins[[paste0("max_", cols[i])]] <- as.numeric(agg[[i]][k])
  bins <- bins[order(bins$chr, bins$bin_start_bp), ]
  rownames(bins) <- NULL
  list(points = points, bins = bins, keep_p = keep_p, bin_bp = bin_bp)
}

# Greedy distance clumping: take the most significant variant, drop every variant within
# window_bp of it on the same chromosome, and repeat while variants pass p_threshold.
#   dt: data.frame with chr, bp and the -log10 p column col.
# Value: the lead rows, most significant first.
distance_clump <- function(dt, col, p_threshold = 5e-8, window_bp = 1e6) {
  dt <- as.data.frame(dt)
  thr <- -log10(p_threshold)
  cand <- which(is.finite(dt[[col]]) & dt[[col]] >= thr)
  if (!length(cand)) return(dt[0, , drop = FALSE])
  cand <- cand[order(-dt[[col]][cand])]
  taken <- logical(length(cand)); alive <- rep(TRUE, length(cand))
  chr <- dt$chr[cand]; bp <- dt$bp[cand]
  for (i in seq_along(cand)) {
    if (!alive[i]) next
    taken[i] <- TRUE
    near <- chr == chr[i] & abs(bp - bp[i]) <= window_bp
    alive[near] <- FALSE
  }
  out <- dt[cand[taken], , drop = FALSE]
  rownames(out) <- NULL
  out
}

# Pooled Stage 1 covariance Sigma0 from a Stage 1 output file (scan.sigma0.tsv.gz: one row
# per quantile level, columns s_<tau>), reordered to the quantile levels taus and symmetrised.
gw_read_sigma0 <- function(path, taus) {
  s0 <- data.table::fread(cmd = paste("gzip -dc", shQuote(path)), showProgress = FALSE)
  scols <- grep("^s_", names(s0), value = TRUE)
  ri <- match(round(taus, 6), round(as.numeric(s0$tau), 6))
  ci <- match(round(taus, 6), round(as.numeric(sub("^s_", "", scols)), 6))
  if (anyNA(ri) || anyNA(ci)) stop("Sigma_0 taus do not align: ", path)
  S <- as.matrix(s0[, scols, with = FALSE])[ri, ci, drop = FALSE]
  (S + t(S)) / 2
}

# Per-chromosome pooled Sigma0 of one Stage 1 GWAS, as exported with the lead profiles
# (results/01_data/lead_profiles/sigma0_<source>.csv).
#   source: Stage 1 GWAS name (the stage1_source column of lead_profiles.csv); u: quantile
#   levels; results_root: the results directory.
# Value: named list of n_tau x n_tau matrices, "1".."22" per chromosome and "mean" (their
#   average, used as the single whitening frame for a trait).
p05_sigma0 <- function(source, u, results_root = file.path(repo_root(), "results")) {
  f <- file.path(results_root, "01_data", "lead_profiles", paste0("sigma0_", source, ".csv"))
  s <- data.table::fread(f, showProgress = FALSE)
  sc <- grep("^s_", names(s), value = TRUE)
  st <- as.numeric(sub("^s_", "", sc))
  ci <- match(round(u, 3), round(st, 3))
  if (anyNA(ci)) stop("sigma0 columns do not cover the tau grid: ", f)
  out <- lapply(split(s, s$chr), function(x) {
    ri <- match(round(u, 3), round(as.numeric(x$tau), 3))
    if (anyNA(ri)) stop("sigma0 rows do not cover the tau grid: ", f)
    M <- as.matrix(x[, sc, with = FALSE])[ri, ci, drop = FALSE]
    (M + t(M)) / 2
  })
  names(out) <- sub("^0+", "", as.character(names(out)))
  out
}

# One lead set from lead_profiles.csv as a profile matrix.
#   P: data.table of lead_profiles.csv; set: lead_set name; u: quantile levels.
# Value: list(meta = per-lead columns (snp_id, chr, bp, alleles, d_resid, ...),
#   B = leads x n_tau matrix of quantile effects beta_tau_<u>, source = Stage 1 GWAS name).
lead_set_matrix <- function(P, set, u) {
  x <- P[P$lead_set == set, ]
  if (!nrow(x)) stop("lead set not in the lead profiles: ", set)
  bc <- paste0("beta_tau_", sprintf("%.3f", u))
  if (!all(bc %in% names(x))) stop("lead profile columns do not match the tau grid")
  B <- as.matrix(x[, bc, with = FALSE]); storage.mode(B) <- "double"
  list(meta = as.data.frame(x[, c("lead_set", "trait", "stage1_source", "snp_id", "chr", "bp",
                                   "effect_allele", "other_allele", "aligned_flip", "d_resid"),
                               with = FALSE]),
       B = B, source = unique(x$stage1_source))
}
