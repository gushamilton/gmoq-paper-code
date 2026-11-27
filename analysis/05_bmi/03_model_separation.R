#!/usr/bin/env Rscript
# Simulation-calibrated comparison of equal-rank BMI models.
# Q_M = d beta' (I-H_M)' Sigma0^-1 (I-H_M) beta with an unweighted OLS
# projection H_M. For each lead, calibrate Q_B-Q_A under both pure-model
# Gaussian truths, matched on the observed truncated mean effect. Common
# random numbers are shared across leads and model pairs within a chromosome.
# pA_obs is the lower-tail probability under A; pB_obs the upper tail under B.
# At nominal 0.05, rejection of additive on log BMI is 570/639; the reverse is 1/639.
# Inputs are the released lead summaries, pooled covariance and model weights.
# The saved results use 20,000 simulations, seed 20260928. No participant data.
# Run: Rscript analysis/05_bmi/03_model_separation.R [--nsim 20000]
# A 2,000-draw run is a smaller smoke check, not the published calibration.
# At 20,000 draws, each simulation matrix has 900,000 doubles (7.2 MB);
# several matrices are live at once. The aggregate-only run needs less than 1 GB.
args <- commandArgs(trailingOnly = TRUE)
argv <- function(k, default = NULL) { i <- which(args == k); if (length(i)) args[i + 1] else default }
repo <- Sys.getenv("GMOQ_REPO", "")
if (!nzchar(repo)) repo <- system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE)
source(file.path(repo, "R", "load.R"))
cfg <- read_config()
OUT <- "05_bmi/model_separation"
NSIM <- as.integer(argv("--nsim", "20000")); SEED <- 20260928L
u <- tau_grid(cfg); nt <- length(u); cw <- trap_weights(u)
P5 <- as.data.table(read_result("01_data/lead_profiles", "lead_profiles"))
WMt <- as.data.table(read_result("05_bmi/model_scores", "weight_matrices"))
SCt <- as.data.table(read_result("05_bmi/model_scores", "lead_scores"))

PRIMARY <- c("additive", "log_additive", "shifted", "both", "logscale_vqtl", "shifted_lognormal_full", "three_components")
wcols <- c("additive", "multiplicative_c0", "multiplicative", "meanlog", "sdlog", "level", "tilt", "bend")
get_W <- function(sc) {
  x <- WMt[trait == "bmi" & scale == sc]
  setNames(lapply(PRIMARY, function(m) {
    y <- x[model == m][order(as.numeric(tau))]
    if (nrow(y) != nt) stop("weight matrix rows: ", sc, " ", m)
    cc <- wcols[vapply(wcols, function(k) all(!is.na(y[[k]])), TRUE)]
    W <- as.matrix(y[, cc, with = FALSE]); storage.mode(W) <- "double"; W
  }), PRIMARY)
}
Kmat <- function(X, Si) { M <- diag(nt) - X %*% solve(crossprod(X), t(X)); t(M) %*% Si %*% M }
qf <- function(B, K) rowSums((B %*% K) * B)
PAIRS <- list(c("additive", "log_additive"), c("additive", "shifted"), c("log_additive", "shifted"))
set.seed(SEED); E <- matrix(rnorm(NSIM * nt), NSIM, nt)

lead_rows <- list(); qall <- list()
for (sc in c("raw", "log")) {
  set_name <- if (sc == "raw") "bmi_A_regenie" else "log_bmi_A_regenie"
  ls <- lead_set_matrix(P5, set_name, u); S0 <- p05_sigma0(ls$source, u)
  meta <- as.data.table(ls$meta); B <- ls$B; W <- get_W(sc)
  stopifnot(nrow(meta) == 639L)
  meta[, `:=`(delta_mu = as.numeric(B %*% cw))]
  dz <- P5[lead_set == set_name]$delta_mu_z[match(meta$snp_id, P5[lead_set == set_name]$snp_id)]
  meta[, delta_mu_z := dz]
  for (k in unique(as.character(meta$chr))) {
    i <- which(as.character(meta$chr) == k); d <- meta$d_resid[i]
    Sk <- S0[[k]]; Si <- solve(Sk); Si <- (Si + t(Si)) / 2
    R <- chol(Sk); ER <- E %*% R                                   # rows ~ N(0, Sigma0)
    Sih <- { e <- eigen(Sk, symmetric = TRUE); e$vectors %*% diag(1 / sqrt(e$values)) %*% t(e$vectors) }
    K <- lapply(W, Kmat, Si = Si)
    Q <- sapply(PRIMARY, function(m) d * qf(B[i, , drop = FALSE], K[[m]]))
    if (length(i) == 1L) Q <- matrix(Q, 1, dimnames = list(NULL, PRIMARY))
    qall[[paste(sc, k)]] <- data.table(scale = sc, snp_id = meta$snp_id[i], delta_mu_z = meta$delta_mu_z[i], Q)
    for (pr in PAIRS) {
      A <- pr[1]; Bm <- pr[2]; wA <- as.numeric(W[[A]]); wB <- as.numeric(W[[Bm]])
      D0 <- K[[Bm]] - K[[A]]; ERD <- ER %*% D0; C <- rowSums(ERD * ER)
      # Targeted GLS angle in this lead's whitened metric (d cancels)
      yA <- Sih %*% wA; yB <- Sih %*% wB
      sinA <- sqrt(max(0, 1 - sum(yA * yB)^2 / (sum(yA^2) * sum(yB^2))))
      for (j in seq_along(i)) {
        ii <- i[j]; dm <- meta$delta_mu[ii]; zdm <- abs(meta$delta_mu_z[ii])
        muA <- wA * dm / sum(cw * wA); muB <- wB * dm / sum(cw * wB)
        aA <- d[j] * sum(muA * (D0 %*% muA)); aB <- d[j] * sum(muB * (D0 %*% muB))
        bA <- sqrt(d[j]) * as.numeric(ERD %*% muA); bB <- sqrt(d[j]) * as.numeric(ERD %*% muB)
        dq <- function(s, a, b) s^2 * a + 2 * s * b + C
        pw <- function(s, dir) {                               # dir A: decide A (dQ large)
          qa <- dq(s, aA, bA); qb <- dq(s, aB, bB)
          if (dir == "A") mean(qa > quantile(qb, 0.95, names = FALSE)) else mean(qb < quantile(qa, 0.05, names = FALSE))
        }
        s80 <- function(dir) {
          if (pw(100, dir) < 0.8) return(Inf)
          lo <- log(1e-3); hi <- log(100)
          for (it in 1:30) { mid <- (lo + hi) / 2; if (pw(exp(mid), dir) >= 0.8) hi <- mid else lo <- mid }
          exp(hi)
        }
        qa1 <- dq(1, aA, bA); qb1 <- dq(1, aB, bB); obs <- Q[j, Bm] - Q[j, A]
        zw <- sqrt(d[j]) * sqrt(sum((Sih %*% muA)^2))            # whitened norm of the pure-A truth
        zwB <- sqrt(d[j]) * sqrt(sum((Sih %*% muB)^2))
        lead_rows[[length(lead_rows) + 1L]] <- data.table(scale = sc, snp_id = meta$snp_id[ii], chr = meta$chr[ii],
          delta_mu_z = meta$delta_mu_z[ii], model_a = A, model_b = Bm,
          lambda_truthA = d[j] * sum(muA * (K[[Bm]] %*% muA)), lambda_truthB = d[j] * sum(muB * (K[[A]] %*% muB)),
          mean_dQ_truthA = mean(qa1), mean_dQ_truthB = mean(qb1), tr_DSigma = mean(C),
          p_correct_truthA = mean(qa1 > 0), p_correct_truthB = mean(qb1 < 0),
          p_nontie_correct_truthA = mean(qa1 > 2), p_nontie_correct_truthB = mean(qb1 < -2),
          power_test_truthA = pw(1, "A"), power_test_truthB = pw(1, "B"),
          z80_test_truthA = s80("A") * zdm, z80_test_truthB = s80("B") * zdm,
          whitened_norm_truthA = zw, whitened_norm_truthB = zwB, angle_whitened_deg = asin(min(1, sinA)) * 180 / pi,
          z80_P10_whitened = 2.8016 / sinA, z80_P10_dmz_truthA = 2.8016 / sinA * zdm / zw,
          z80_P10_dmz_truthB = 2.8016 / sinA * zdm / zwB,
          dQ_obs = obs, pA_obs = mean(qa1 <= obs), pB_obs = mean(qb1 >= obs))
      }
    }
  }
}
LD <- rbindlist(lead_rows); QA <- rbindlist(qall)
# cross-check Q against saved lead scores
x <- SCt[set == "primary", .(scale, snp_id, model, Q_full)]
qm <- melt(QA, id.vars = c("scale", "snp_id", "delta_mu_z"), variable.name = "model", value.name = "Q_here")
qm[, model := as.character(model)]
XC <- merge(qm, x, by = c("scale", "snp_id", "model"))
XCS <- XC[, .(n = .N, max_rel_diff = max(abs(Q_here - Q_full) / pmax(Q_full, 1e-12))), by = .(scale, model)]
print(XCS)
if (nrow(XC) != nrow(qm) || max(XCS$max_rel_diff) > 1e-6) stop("Q does not reproduce saved lead scores")

bins <- c(0, 7, 10, 20, 40, Inf)
LD[, strength_bin := cut(abs(delta_mu_z), bins, include.lowest = TRUE)]
summ <- function(x) x[, .(n_leads = .N,
  median_dmz = median(abs(delta_mu_z)),
  frac_power80_test_truthA = mean(power_test_truthA >= 0.8), frac_power80_test_truthB = mean(power_test_truthB >= 0.8),
  frac_pcorrect80_truthA = mean(p_correct_truthA >= 0.8), frac_pcorrect80_truthB = mean(p_correct_truthB >= 0.8),
  frac_pnontie80_truthA = mean(p_nontie_correct_truthA >= 0.8), frac_pnontie80_truthB = mean(p_nontie_correct_truthB >= 0.8),
  median_power_test_truthA = median(power_test_truthA), median_power_test_truthB = median(power_test_truthB),
  median_z80_test_truthA = median(z80_test_truthA), median_z80_test_truthB = median(z80_test_truthB),
  median_z80_P10_dmz_truthA = median(z80_P10_dmz_truthA), median_z80_P10_whitened = median(z80_P10_whitened),
  median_angle_whitened_deg = median(angle_whitened_deg),
  median_ratio_whitened_to_dmz = median(whitened_norm_truthA / abs(delta_mu_z)),
  median_tr_DSigma = median(tr_DSigma),
  frac_obs_rejects_B = mean(pB_obs < 0.05), frac_obs_rejects_A = mean(pA_obs < 0.05),
  frac_obs_rejects_neither = mean(pA_obs >= 0.05 & pB_obs >= 0.05), frac_obs_rejects_both = mean(pA_obs < 0.05 & pB_obs < 0.05),
  frac_obs_absdQ_ge2 = mean(abs(dQ_obs) >= 2), median_obs_dQ = median(dQ_obs),
  q10_obs_dQ = quantile(dQ_obs, 0.1), q90_obs_dQ = quantile(dQ_obs, 0.9)), by = c("scale", "model_a", "model_b", intersect(names(x), "strength_bin"))]
PS <- summ(LD[, !"strength_bin"]); PB <- summ(LD); setorder(PB, scale, model_a, model_b, strength_bin)
print(PS)


write_result(LD, OUT, "lead_pair_power")
write_result(PS, OUT, "pair_power_summary")
write_result(PB, OUT, "pair_power_by_strength")
write_result(XCS, OUT, "q_crosscheck")
