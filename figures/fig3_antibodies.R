#!/usr/bin/env Rscript
# Figure 3: decomposition of antibody associations into membership and level effects
#
# What: antibody levels against JC and Merkel cell polyomavirus (JCV, MCV) VP1 are fitted as a
#   two-component skew-normal mixture. Stage 2 then splits each lead variant's quantile
#   profile into effects on the lower- and upper-component means and on membership of the
#   upper component, set beside REGENIE on the mean, on seropositivity and on IgG among
#   seropositives.
#   a: phenotype histograms with the fitted mixtures.
#   b-d: three lead variants (JCV HLA, JCV FUT2, MCV HLA): estimates and Stage 1 profiles.
# Paper: Fig. 3
# Input: results/04_antibodies/decomposition/{phenotype_histograms,mixture_density,
#   mixture_components,lead_decomposition,example_profiles}.csv,
#   results/04_antibodies/serostatus_gwas/serostatus_leads.csv
# Output: paper_assets/fig3_antibodies.{png,pdf}, _source.csv, _caption.md
# Runs: locally from results/ in about 3 seconds

suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(data.table)})
source(file.path(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE), "figures", "theme.R"))

cfg <- read_config()
trait_lab <- c(mcv_vp1_log_all = "MCV VP1", jcv_vp1_log_all = "JCV VP1")
# Each antibody is reported with its weight design from config.yaml ('model' for JCV VP1,
# 'hybrid' for the bimodal MCV VP1).
antibody <- function(name) pick_design(read_result("04_antibodies/decomposition", name), "trait", cfg = cfg)
col_lower <- "#4a7fb5"; col_upper <- "#c0392b"
part_col <- c(total = "grey10", lower = col_lower, upper = col_upper, gamma = paper_palette[["mixture"]])

## 1. Panel a: distributions and fitted mixtures ---------------------------------------------
hist <- antibody("phenotype_histograms")[is.finite(bin_lower) & is.finite(bin_upper) & !is.na(density)]
mix_density <- antibody("mixture_density")
for (d in list(hist, mix_density)) d[, trait_f := factor(trait_lab[trait], levels = trait_lab)]
mix_long <- melt(mix_density, id.vars = c("trait_f", "x"), measure.vars = c("density_total", "density_lower", "density_upper"),
                 variable.name = "part", value.name = "density")
mix_long[, part := factor(sub("density_", "", part), levels = c("total", "lower", "upper"))]
mix_long <- mix_long[density > 1e-3 | part == "total"]
pA <- ggplot() +
  geom_rect(data = hist, aes(xmin = bin_lower, xmax = bin_upper, ymin = 0, ymax = density), fill = "grey84", colour = NA) +
  geom_line(data = mix_long[part != "total"], aes(x, density, colour = part), linewidth = 0.6, lineend = "round") +
  geom_line(data = mix_long[part == "total"], aes(x, density, colour = part), linewidth = 0.85, lineend = "round") +
  facet_wrap(~trait_f, nrow = 1) +
  scale_colour_manual(values = part_col[c("total", "lower", "upper")], labels = c(total = "mixture", lower = "lower", upper = "upper"),
                      breaks = c("total", "lower", "upper"), name = NULL) +
  scale_x_continuous(limits = c(1, 10.2), breaks = seq(2, 10, 2)) +
  labs(x = "log antibody level (MFI)", y = "Density", tag = "a") +
  theme_paper() + legend_inside(0.995, 0.97, c(1, 1))

## 2. Panels b-d: three lead variants -----------------------------------------------------
examples <- data.table(trait = c("jcv_vp1_log_all", "jcv_vp1_log_all", "mcv_vp1_log_all"),
                       snp_id = c("rs9271055", "rs681343", "rs9271055"), gene = c("HLA", "FUT2", "HLA"))
leads <- antibody("lead_decomposition"); profiles <- antibody("example_profiles")
serostatus <- as.data.table(read_result("04_antibodies/serostatus_gwas", "serostatus_leads"))
forest_col <- c(regenie = paper_palette[["regenie"]], mean = paper_palette[["gmoq"]], lower = col_lower, upper = col_upper,
                gamma = paper_palette[["mixture"]])
profile_lab <- c(total = "total", lower = "lower", upper = "upper", gamma = "membership")
sero_lab <- c(seropositivity = "REGENIE: seropositivity (log OR)", seropositive_igg = "REGENIE: IgG in seropositives")
rows <- list(); src_ex <- list()
for (k in seq_len(nrow(examples))) {
  x <- leads[trait == examples$trait[k] & snp_id == examples$snp_id[k]]; stopifnot(nrow(x) == 1)
  s <- sign(x$regenie_beta)             # align to the allele with a positive REGENIE mean effect
  sr <- serostatus[trait == examples$trait[k] & snp_id == examples$snp_id[k]]
  fd <- rbindlist(lapply(names(sero_lab), function(a)
    data.table(label = sero_lab[[a]], est = s * sr[analysis == a]$beta, se = sr[analysis == a]$se, ck = "regenie", block = "r")))
  # The mixture's gamma is the log odds of the lower component, so -gamma is membership of
  # the upper component.
  fd <- rbind(fd,
    data.table(label = "REGENIE: mean log MFI (all)", est = s * x$regenie_beta, se = x$regenie_se, ck = "regenie", block = "r"),
    data.table(label = "G-MoQ: mean", est = s * x$mean_effect, se = x$mean_se, ck = "mean", block = "g"),
    data.table(label = "G-MoQ: lower-component mean", est = s * x$beta_left, se = x$beta_left_se, ck = "lower", block = "g"),
    data.table(label = "G-MoQ: upper-component mean", est = s * x$beta_right, se = x$beta_right_se, ck = "upper", block = "g"),
    data.table(label = "G-MoQ: membership (log odds upper)", est = -s * x$gamma, se = x$gamma_se, ck = "gamma", block = "g"))
  pr <- profiles[trait == examples$trait[k] & snp_id == examples$snp_id[k]]
  obs <- pr[, .(tau, beta = s * beta_stage1, se = se_stage1)]
  fits <- rbind(pr[, .(tau, series = "total", value = s * fitted_total)], pr[, .(tau, series = "lower", value = s * fitted_beta_left)],
                pr[, .(tau, series = "upper", value = s * fitted_beta_right)], pr[, .(tau, series = "gamma", value = s * fitted_gamma)])
  pf <- snp_forest(fd, forest_col, xlab = "Per-allele effect (95% CI)")
  pp <- snp_profile(obs, fits, part_col, labels = profile_lab, legend_pos = c(0.02, 0.99), legend_just = c(0, 1),
                    ylab = "Stage 1 beta_tau (95% CI)", headroom = if (k == 1) 0.5 else 0)
  if (k > 1) pp <- pp + theme(legend.position = "none")
  header <- sprintf("%s, %s %s", trait_lab[[examples$trait[k]]], examples$gene[k], examples$snp_id[k])
  rows[[k]] <- (pf + labs(tag = letters[1 + k]) + ggtitle(header) +
                  theme(plot.title = element_text(size = 7, face = "plain", hjust = 0,
                                                  margin = margin(b = 2))) | pp) +
    plot_layout(widths = c(1, 1.35))
  src_ex[[k]] <- rbind(fd[, .(panel = LETTERS[k + 1], trait = examples$trait[k], group = examples$snp_id[k], series = label, x = NA_real_, y = est, se)],
                       obs[, .(panel = LETTERS[k + 1], trait = examples$trait[k], group = examples$snp_id[k], series = "stage1", x = tau, y = beta, se)],
                       fits[, .(panel = LETTERS[k + 1], trait = examples$trait[k], group = examples$snp_id[k], series = paste0("fit_", series), x = tau, y = value, se = NA_real_)])
}

## 3. Assemble and save ---------------------------------------------------------------------
fig <- pA / rows[[1]] / rows[[2]] / rows[[3]] + plot_layout(heights = c(0.85, 1, 1, 1))
src <- rbind(
  hist[, .(panel = "A", trait, group = "histogram", series = "density", x = (bin_lower + bin_upper) / 2, y = density, se = NA_real_)],
  mix_long[, .(panel = "A", trait = as.character(trait_f), group = "mixture", series = as.character(part), x, y = density, se = NA_real_)],
  rbindlist(src_ex))
cap <- paste(
  "Figure 3. Decomposition of antibody associations into membership and level effects.",
  "UK Biobank JC and Merkel cell polyomavirus VP1 IgG, log mean fluorescence intensity; n = 9,145.",
  "**a**, Phenotype histograms with the fitted two-component skew-normal mixture that defines the Stage 2 weights:",
  "total density (mixture) and the lower and upper components. MCV VP1 is bimodal, so its components can be read as",
  "seronegative-like and seropositive-like groups; JCV VP1 is unimodal and its components describe shape rather than",
  "exposure groups. **b**-**d**, Three lead variants: JCV VP1 at HLA (rs9271055), JCV VP1 at FUT2 (rs681343) and MCV",
  "VP1 at HLA (rs9271055). Left, per-allele effects with 95% confidence intervals from REGENIE on all-sample log MFI,",
  "on assay-defined seropositivity (log odds) and on IgG among seropositive participants, beside the G-MoQ mean and",
  "the three mixture parameters: lower-component mean, upper-component mean and the log odds of membership of the",
  "upper component. Right, Stage 1 quantile effects (grey, 95% confidence intervals) with the fitted profile (black)",
  "and each parameter's contribution. At JCV, HLA acts mainly through membership, and at MCV it also moves both component",
  "means; FUT2 shifts both component means with little",
  "effect on membership. Stage 2 is ordinary least squares with the pooled Stage 1 covariance (Methods).",
  "MFI, mean fluorescence intensity.")
save_figure(fig, "fig3_antibodies", src, cap, width = 183, height = 230,
            results_used = c("04_antibodies/decomposition", "04_antibodies/serostatus_gwas"))
