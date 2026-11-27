#!/usr/bin/env Rscript
# Fig 2: testosterone recovery of sex-specific effects from a sex-blind quantile GWAS.
# Layout after the author's review (28 Sep 2026), old Fig 2 style:
#   (A) top: G-MoQ component vs sex-stratified reference REGENIE at clumped reference lead SNPs;
#       compute their agreement directly, keeping all-significant-SNP statistics separate.
#   (B) matched LDSC rg by component and shrinkage (old Fig 2E tile style); (C) cost per additional SNP;
#   (D, E) two exemplar SNPs picked by the pre-written selection rule (P08 exemplar_selected, unshrunk,
#   hybrid): the mixed-component SNP (rule category shared_concordant: female, male and gamma all
#   |z| > 1.96) and the male-specific SNP. Each exemplar: LEFT forest of sex-stratified reference vs
#   G-MoQ per-allele estimates (female, male, mixing), RIGHT Stage 1 profile + fitted curves.
# Use the reported hybrid design in every arm, matching Supplementary Table S3.
# Retain the dashed outline marking the unshrunk near-empty-gap example.
# Reads results/P08_testosterone (+ P08b_exemplar_sex_reference when present), P07 weight_density and
# aux_speed_profile per_snp_costs (panel C).
suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(data.table)})
.repo <- system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE)
source(file.path(.repo, "figures", "theme.R"))

cfg <- read_config()
arm_lab <- c(testosterone = "No shrinkage", testosterone_shrink_0_2 = "Shrink 0.2",
             testosterone_shrink_0_5 = "Shrink 0.5", testosterone_shrink_0_7 = "Shrink 0.7")
comp_col <- c(female = paper_palette[["female"]], male = paper_palette[["male"]],
              gamma = paper_palette[["mixture"]], total = "grey10")
D0 <- mixture_design(cfg, "testosterone")
rr <- function(name, all = FALSE) { x <- as.data.table(read_result("03_testosterone/recovery", name))
  if (all) x else if ("arm" %in% names(x)) pick_design(x, "arm", cfg = cfg) else pick_design(x, fit = "testosterone", cfg = cfg) }
flag_arm <- "testosterone"   # unshrunk near-empty-gap example

# ---- A: lead-level G-MoQ vs reference (top) --------------------------------------------------
cl <- rr("concordance_leads")[reference == "combined" & component %in% c("female", "male")]
cl[, arm_f := factor(arm_lab[arm], levels = arm_lab)]
cc <- rr("concordance")
ccl <- cl[, .(n = .N, pearson_beta = cor(beta_regenie, beta_gmoq),
              slope_gmoq_on_regenie = unname(coef(lm(beta_gmoq ~ beta_regenie))[2])),
          by = .(arm, component)]
ccl[, arm_f := factor(arm_lab[arm], levels = arm_lab)]
cl[, comp_f := factor(component, levels = c("female", "male"), labels = c("Female effect", "Male effect"))]
ccl[, comp_f := factor(component, levels = c("female", "male"), labels = c("Female effect", "Male effect"))]
flag <- CJ(arm_f = factor(arm_lab[[flag_arm]], levels = arm_lab), comp_f = levels(cl$comp_f))
flag[, comp_f := factor(comp_f, levels = levels(cl$comp_f))]
flag[, `:=`(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)]
pA <- ggplot(cl, aes(beta_regenie, beta_gmoq, colour = component)) +
  geom_rect(data = flag, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), inherit.aes = FALSE,
            fill = NA, colour = "grey45", linetype = "22", linewidth = 0.35) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.3, linetype = 2, colour = "grey45") +
  geom_point(size = 0.7, alpha = 0.65, stroke = 0) +
  facet_grid(comp_f ~ arm_f) +
  scale_x_continuous(breaks = c(-0.05, 0, 0.05)) +
  scale_colour_manual(values = comp_col, guide = "none") +
  labs(x = "Sex-stratified reference REGENIE beta", y = "G-MoQ Stage 2 beta") +
  theme_paper() + theme(strip.text.y.right = element_text(angle = -90, hjust = 0.5), panel.spacing = unit(4, "pt"))

# ---- B: LDSC rg matrix (old Fig 2E style) -----------------------------------------------------
rg <- rr("ldsc_rg")
rg[, arm_f := factor(arm_lab[arm], levels = arm_lab)]
rg[, comp_f := factor(component, levels = rev(c("female", "male", "sex")),
                      labels = rev(c("Female effect", "Male effect", "Upper membership\nvs sex GWAS")))]
pB <- ggplot(rg, aes(arm_f, comp_f, fill = rg)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey10"), guide = "none") +
  scale_fill_gradient(low = "#f1e6d8", high = "#2f6690", limits = c(min(0.7, min(rg$rg)), 1),
                      oob = scales::squish, name = "rg") +  # rg can exceed 1 by sampling error (1.0003)
  labs(x = NULL, y = NULL) +
  theme_paper() + theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(),
                        legend.position = "right", legend.key.height = unit(1.2, "lines"),
                        legend.key.width = unit(0.5, "lines"))

# ---- C: cost per additional SNP (speed_rerun/results/component_costs.csv) ---------------------
# One core, warm session; one-off per-trait costs (baseline mixture fit, Stage 1 density) excluded.
pc <- as.data.table(read_result("08_simulations/speed_benchmark", "component_costs"))
c_lv <- c("Direct mixture ML", "G-MoQ Stage 1 + Stage 2", "Stage 2 alone (given Stage 1)")
c_col <- setNames(c(paper_palette[["mixture"]], paper_palette[["gmoq"]], "#7fb3d5"), c_lv)
tm <- rbind(pc[, .(n, method = c_lv[1], s = direct_seconds_per_variant)],
            pc[, .(n, method = c_lv[2], s = stage1_seconds_per_extra_variant + stage2_seconds_per_variant_given_summaries)],
            pc[, .(n, method = c_lv[3], s = stage2_seconds_per_variant_given_summaries)])
tm[, method := factor(method, levels = c_lv)]
c400 <- tm[n == 4e5][order(method)]
pC <- ggplot(tm, aes(n, s, colour = method)) +
  geom_line(linewidth = 0.7, lineend = "round") + geom_point(size = 1.3) +
  scale_x_log10(breaks = c(1e4, 5e4, 4e5), labels = c("10k", "50k", "400k"), expand = expansion(mult = c(0.05, 0.28))) +
  scale_y_log10(breaks = 10^(-5:1), labels = c("10 us", "0.1 ms", "1 ms", "10 ms", "0.1 s", "1 s", "10 s"), limits = c(3e-6, 25)) +
  scale_colour_manual(values = c_col, breaks = c_lv, labels = c("Mixture ML", "Stages 1 + 2", "Stage 2"), name = NULL) +
  labs(x = "Sample size", y = "Seconds per SNP") +
  theme_paper() + legend_inside(0.02, 0.66, c(0, 1))

# ---- D, E: exemplars (single-SNP layout) ------------------------------------------------------
es <- rr("exemplar_selected")
ee <- rr("exemplar_estimates")[arm == flag_arm]
ep <- rr("exemplar_profiles")[arm == flag_arm]
sexref <- if (file.exists(file.path(.repo, "results", "03_testosterone", "exemplars", "exemplar_sex_reference.csv")))
  as.data.table(read_result("03_testosterone/exemplars", "exemplar_sex_reference")) else NULL
pick <- c(shared_concordant = "Mixed-component signal", male_specific = "Male-specific signal")
fcols <- c(ref = "grey25", female = comp_col[["female"]], male = comp_col[["male"]], gamma = comp_col[["gamma"]])
pcols <- c(total = "grey10", female = comp_col[["female"]], male = comp_col[["male"]], gamma = comp_col[["gamma"]])
plabs <- c(total = "total", female = "female", male = "male", gamma = "membership")
ex_rows <- list(); src_ex <- list()
for (k in seq_along(pick)) {
  cat_k <- names(pick)[k]; s <- es[category == cat_k]
  e <- ee[snp_id == s$snp_id]; pr <- ep[snp_id == s$snp_id]
  sx <- if (!is.null(sexref)) sexref[snp_id == s$snp_id & weight_design == D0] else NULL
  fd <- rbind(
    data.table(label = "Female: reference REGENIE", est = e$regenie_combined_female_beta, se = e$regenie_combined_female_se, ck = "ref", block = "f"),
    data.table(label = "Female: G-MoQ", est = e$beta_left, se = e$beta_left_se, ck = "female", block = "f"),
    data.table(label = "Male: reference REGENIE", est = e$regenie_combined_male_beta, se = e$regenie_combined_male_se, ck = "ref", block = "m"),
    data.table(label = "Male: G-MoQ", est = e$beta_right, se = e$beta_right_se, ck = "male", block = "m"),
    if (!is.null(sx) && nrow(sx)) data.table(label = "Membership: sex GWAS", est = sx$sex_beta, se = sx$sex_se, ck = "ref", block = "g"),
    # Stored gamma is lower-component log odds; display manuscript upper-component log odds.
    data.table(label = "Membership: G-MoQ", est = -e$gamma, se = e$gamma_se, ck = "gamma", block = "g"))
  obs <- pr[, .(tau, beta = beta_stage1, se = se_stage1)]
  fits <- rbind(pr[, .(tau, series = "total", value = fitted_total)], pr[, .(tau, series = "female", value = fitted_beta_left)],
                pr[, .(tau, series = "male", value = fitted_beta_right)], pr[, .(tau, series = "gamma", value = fitted_gamma)])
  # display range from the taus away from the empty gap at tau = 0.5 (the unshrunk Stage 1 SE there is huge)
  away <- pr[abs(tau - 0.5) > 0.03]
  yl <- range(c(away$beta_stage1 - 1.96 * away$se_stage1, away$beta_stage1 + 1.96 * away$se_stage1, pr$fitted_total[abs(pr$tau - 0.5) > 0.03]))
  yl <- yl + c(-0.08, 0.08) * diff(yl)
  pf <- snp_forest(fd, fcols, xlab = "Per-allele estimate (95% CI)")
  pp <- snp_profile(obs, fits, pcols, labels = plabs, ylim = yl,
                    legend_pos = if (k == 1) c(0.02, 0.98) else c(0.02, 0.03), legend_just = if (k == 1) c(0, 1) else c(0, 0)) +
    labs(y = "Stage 1 beta_tau (95% CI)")
  ex_rows[[k]] <- (pf + labs(tag = letters[3 + k]) | pp) +
    plot_layout(widths = c(1, 1.45))
  src_ex[[k]] <- rbind(fd[, .(panel = LETTERS[3 + k], group = s$snp_id, series = label, x = NA_real_, y = est, se)],
                       obs[, .(panel = LETTERS[3 + k], group = s$snp_id, series = "stage1", x = tau, y = beta, se)],
                       fits[, .(panel = LETTERS[3 + k], group = s$snp_id, series = paste0("fit_", series), x = tau, y = value, se = NA_real_)])
}
pA <- pA + labs(tag = "a"); pB <- pB + labs(tag = "b"); pC <- pC + labs(tag = "c")
p <- pA / ((pB | pC) + plot_layout(widths = c(1.15, 1))) / ex_rows[[1]] / ex_rows[[2]] +
  plot_layout(heights = c(1.25, 0.8, 1, 1))

src <- rbind(
  cl[, .(panel = "A", group = paste(arm, component), series = snp_id, x = beta_regenie, y = beta_gmoq, se = se_gmoq)],
  rg[, .(panel = "B", group = arm, series = component, x = NA_real_, y = rg, se = rg_se)],
  tm[, .(panel = "C", group = as.character(method), series = "seconds_per_snp", x = n, y = s, se = NA_real_)],
  rbindlist(src_ex))
u0 <- cc[arm == "testosterone" & reference == "combined"]
nall <- u0[component == "female" & snp_set == "all", n][1]
g0 <- function(comp, set, col) u0[component == comp & snp_set == set][[col]][1]
wd0 <- as.data.table(read_result("02_mixtures/stage2_mixtures", "weight_density"))[fit == "testosterone"]
wd0 <- wd0[which.min(abs(tau - 0.5))]
sx_note <- if (is.null(sexref)) paste(
  "The chromosomal-sex GWAS reference for the mixing parameter at the exemplar SNPs is not yet exported",
  "(pipeline step P08b, pending), so the mixing row shows G-MoQ only.") else
  sprintf("The mixing row compares -gamma with the chromosomal-sex GWAS on the same %s participants (scales differ).",
          format(cfg$p08$sex_reference_n, big.mark = ","))
cap <- paste(
  "Figure 2. Recovery of sex-specific testosterone effects from a sex-blind analysis.",
  "n = 399,710; the sex-stratified reference is REGENIE in the same 199,572 female and 200,138 male participants.",
  "Hybrid membership weights are used in all arms (Methods).",
  "**a**, G-MoQ female and male estimates against the sex-stratified reference at clumped reference lead variants",
  "(two-sided P < 5 x 10-8), by degree of shrinkage; the dashed line is y = x.",
  sprintf("For no shrinkage, shrinkage 0.2, 0.5 and 0.7, respectively, female lead counts are %s, Pearson correlations are %s and ordinary least-squares slopes are %s; male lead counts are %s, correlations are %s and slopes are %s.",
          paste(ccl[component == "female"]$n, collapse = ", "),
          paste(sprintf("%.3f", ccl[component == "female"]$pearson_beta), collapse = ", "),
          paste(sprintf("%.2f", ccl[component == "female"]$slope_gmoq_on_regenie), collapse = ", "),
          paste(ccl[component == "male"]$n, collapse = ", "),
          paste(sprintf("%.3f", ccl[component == "male"]$pearson_beta), collapse = ", "),
          paste(sprintf("%.2f", ccl[component == "male"]$slope_gmoq_on_regenie), collapse = ", ")),
  "The all-significant-SNP statistics in Supplementary Table S3 use 1,578 female and 3,048 male reference-significant SNPs before clumping.",
  "**b**, LD score regression genetic correlation of each G-MoQ component with its",
  "reference: sex-stratified GWAS for the female and male components, and a GWAS of chromosomal sex in the same",
  "participants for membership. **c**, Computational cost per additional variant on one core against sample size, for",
  "six-start direct mixture maximum likelihood, G-MoQ Stages 1 and 2, and Stage 2 alone given the Stage 1 output.",
  sprintf(paste("At 400,000 participants these costs are %.3f s, %.3f ms and %.2f microseconds per variant, respectively.",
                "Stage 1 marginal cost is the difference between 300- and 100-variant panel times divided by 200;",
                "Stage 2 is measured on 100,000 columns tiled from 100 simulated variants; direct ML averages five variants.",
                "One-off per-trait costs and storage input/output are excluded."), c400$s[1], c400$s[2] * 1000, c400$s[3] * 1e6),
  "**d**,**e**, Two exemplar loci selected in advance: rs1260326, with effects on all three",
  "parameters, and rs7097461, with a male-specific effect. Left, per-allele estimates with 95% confidence intervals (location on log testosterone, membership on the log-odds scale);",
  "right, Stage 1 quantile effects (grey, 95% confidence intervals) with the fitted profile (black) and each",
  "parameter's contribution. Few participants lie between the modes, so effects near tau = 0.5 are very imprecise and",
  "values beyond the axis are drawn at its edge as triangles. Membership is oriented to the upper (male) component.")

# Export the three denominators separately; the plotted lead statistics are recalculated
# from exactly the same rows as the scatterplot, rather than relabelled all-SNP estimates.
stats_source <- rbind(
  cc[reference == "combined" & component %in% c("female", "male"),
     .(arm, component, snp_set, weight_design, n, pearson_beta, slope_gmoq_on_regenie)],
  ccl[, .(arm, component, snp_set = "plotted_clumped_reference_leads",
          weight_design = mixture_design(cfg, arm), n, pearson_beta, slope_gmoq_on_regenie)])
fwrite(stats_source, file.path(assets_dir, "fig2_testosterone_recovery_statistics.csv"))
save_figure(p, "fig2_testosterone_recovery", src, cap, width = 183, height = 235,
            results_used = c("03_testosterone/recovery/exemplar_profiles", "03_testosterone/recovery/exemplar_selected",
                             "03_testosterone/recovery/exemplar_estimates", "03_testosterone/recovery/concordance_leads",
                             "03_testosterone/recovery/concordance", "03_testosterone/recovery/ldsc_rg", "08_simulations/speed_benchmark/component_costs",
                             "02_mixtures/stage2_mixtures/weight_density",
                             if (!is.null(sexref)) "03_testosterone/exemplars/exemplar_sex_reference"))
