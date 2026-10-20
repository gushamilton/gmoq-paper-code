#!/usr/bin/env Rscript
# Supplementary Figure S9: genome-wide associations with JCV and MCV VP1 antibody levels
#
# What: the genome-wide view behind Fig. 3. The G-MoQ joint (Hotelling) test of the three
#   mixture parameters is set against REGENIE on log antibody level, as mirrored Manhattan
#   plots, with quantile-quantile plots of the joint test, the G-MoQ mean-effect test and
#   REGENIE. Only binned maxima and variants with joint P < 1e-5 are stored and drawn.
#   a: mirrored Manhattan plots (joint test above the axis, REGENIE below), MCV and JCV VP1.
#   b: quantile-quantile plots per antibody.
# Paper: Supplementary Fig. S9
# Input: results/04_antibodies/decomposition/{manhattan_points_*,manhattan_bins_*,qq_*,
#   qq_lambda,overlap}.csv
# Output: paper_assets/figS09_antibody_manhattan_qq.{png,pdf}, _source.csv, _caption.md
# Runs: locally from results/ in about 2 seconds

suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(data.table)})
source(file.path(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE), "figures", "theme.R"))

traits <- c(mcv_vp1_log_all = "MCV VP1", jcv_vp1_log_all = "JCV VP1")
trait_order <- names(traits)
Y_CAP <- 60                                      # -log10 P drawn at most to 60
SIG <- -log10(5e-8)
# Each antibody is drawn with its reported weight design from config.yaml ('model' for JCV VP1,
# 'hybrid' for MCV VP1).
cfg <- read_config()
design <- setNames(mixture_design(cfg, trait_order), trait_order)
antibody <- function(name) as.data.table(read_result("04_antibodies/decomposition", name))

## 1. Load ----------------------------------------------------------------------------------
dat <- lapply(trait_order, function(trait) list(
  pts = antibody(paste0("manhattan_points_", trait, "_", design[[trait]])),
  bins = antibody(paste0("manhattan_bins_", trait, "_", design[[trait]]))))
names(dat) <- trait_order

## 2. Panel a: mirrored Manhattan plots ---------------------------------------------------
# Cumulative genome position, from the bin table (same variant set for both antibodies), so
# chromosome offsets line up across panels.
chr_off <- dat[[1]]$bins[, .(chr_len = max(bin_start_bp) + 1e6), by = chr][order(chr)]
chr_off[, offset := cumsum(shift(chr_len, fill = 0))]
add_pos <- function(dt, bp_col) {
  dt <- merge(dt, chr_off[, .(chr, offset)], by = "chr")
  dt[, pos := get(bp_col) + offset]
  dt
}
axis_breaks <- chr_off[, .(chr, mid = offset + chr_len / 2)]
cap_val <- function(x) pmin(x, Y_CAP)
is_cap <- function(x) x > Y_CAP
manhattan_data <- function(trait) {
  bins <- add_pos(copy(dat[[trait]]$bins), "bin_start_bp")
  pts <- add_pos(copy(dat[[trait]]$pts), "bp")
  bins[, chr_parity := chr %% 2]
  # Grey background: the largest -log10 P in each 1-Mb bin; joint test up, REGENIE down.
  bg <- rbind(
    bins[, .(pos, y = cap_val(max_mlog10p_joint), capped = is_cap(max_mlog10p_joint), chr_parity, side = "joint")],
    bins[, .(pos, y = -cap_val(max_mlog10p_regenie), capped = is_cap(max_mlog10p_regenie), chr_parity, side = "regenie")])
  # Coloured foreground: variants with joint P < 1e-5, on both sides.
  fg <- rbind(
    pts[, .(pos, y = cap_val(mlog10p_joint), capped = is_cap(mlog10p_joint), side = "joint")],
    pts[, .(pos, y = -cap_val(mlog10p_regenie), capped = is_cap(mlog10p_regenie), side = "regenie")])
  bg[, trait_lab := traits[[trait]]]; fg[, trait_lab := traits[[trait]]]
  list(bg = bg, fg = fg)
}
md <- lapply(trait_order, manhattan_data)
bg <- rbindlist(lapply(md, `[[`, "bg")); fg <- rbindlist(lapply(md, `[[`, "fg"))
bg[, trait_lab := factor(trait_lab, levels = traits[trait_order])]
fg[, trait_lab := factor(trait_lab, levels = traits[trait_order])]
# Colour is redundant with the side of the axis, so there is no key.
pMan <- ggplot() +
  geom_point(data = bg[capped == FALSE], aes(x = pos, y = y, colour = factor(chr_parity)),
             size = 0.3, alpha = 0.55, show.legend = FALSE) +
  geom_hline(yintercept = c(SIG, -SIG), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.3) +
  geom_point(data = fg[capped == FALSE], aes(x = pos, y = y, colour = side), size = 0.7, show.legend = FALSE) +
  geom_point(data = fg[capped == TRUE], aes(x = pos, y = ifelse(y > 0, Y_CAP, -Y_CAP), colour = side),
             shape = 17, size = 1.3, show.legend = FALSE) +
  scale_colour_manual(values = c(`0` = "grey75", `1` = "grey55",
                                 joint = unname(paper_palette["gmoq"]),
                                 regenie = unname(paper_palette["regenie"]))) +
  scale_x_continuous(breaks = axis_breaks$mid, labels = axis_breaks$chr, expand = expansion(mult = 0.01)) +
  scale_y_continuous(limits = c(-Y_CAP, Y_CAP) * 1.05, breaks = seq(-Y_CAP, Y_CAP, by = 20),
                     labels = function(v) abs(v)) +
  facet_wrap(~trait_lab, ncol = 1) +
  labs(x = "Chromosome", y = expression(-log[10]~italic(P)), tag = "a") +
  theme_paper() +
  theme(axis.text.x = element_text(size = 5))

## 3. Panel b: quantile-quantile plots ----------------------------------------------------
qq <- rbindlist(lapply(trait_order, function(tr) antibody(paste0("qq_", design[[tr]]))[trait == tr]))
lambda <- pick_design(antibody("qq_lambda"), "trait", cfg = cfg)
qq[, test_label := factor(c(joint = "Joint", mean = "Mean", regenie = "REGENIE")[test],
                          levels = c("Joint", "Mean", "REGENIE"))]
qq[, trait_lab := factor(traits[trait], levels = traits[trait_order])]
pQQ <- ggplot(qq, aes(expected_mlog10p, observed_mlog10p, colour = test_label)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  geom_point(size = 0.5, alpha = 0.7) +
  facet_wrap(~trait_lab, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = c(Joint = unname(paper_palette["gmoq"]),
                                 Mean = unname(paper_palette["mixture"]),
                                 REGENIE = unname(paper_palette["regenie"])), name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 1.4, alpha = 1))) +
  labs(x = expression(Expected~-log[10]~italic(P)), y = expression(Observed~-log[10]~italic(P)), tag = "b") +
  theme_paper() + legend_inside(0.01, 0.98, c(0, 1)) +
  theme(panel.spacing.x = unit(1.2, "lines"))

## 4. Numbers for the legend ----------------------------------------------------------------
# Loci within 500 kb found by each method; the legend says every locus is shared.
overlap <- pick_design(antibody("overlap"), "trait", cfg = cfg)
ov <- overlap[level == "locus_500kb" & gmoq_test %in% c("gmoq_joint", "gmoq_mean")]
stopifnot(ov[, all(shared == n_gmoq & shared == n_regenie)])
n_loci <- ov[gmoq_test == "gmoq_joint", setNames(n_regenie, trait)]
lambda_txt <- function(tr) {
  l <- lambda[trait == tr]
  sprintf("%s VP1, joint %.3f, mean %.3f, REGENIE %.3f", toupper(sub("_vp1_log_all", "", tr)),
          l$lambda_joint, l$lambda_mean, l$lambda_regenie)
}
n_snps <- unique(lambda$n_snps); stopifnot(length(n_snps) == 1)

## 5. Assemble and save ---------------------------------------------------------------------
fig <- pMan / pQQ + plot_layout(heights = c(2.1, 1))
src <- rbindlist(lapply(trait_order, function(tr) cbind(trait = tr, dat[[tr]]$pts)), fill = TRUE)
cap <- paste(
  "Figure S9. Genome-wide associations with JC and Merkel cell polyomavirus VP1 antibody levels.",
  sprintf("UK Biobank log MFI; n = 9,145 participants and %s variants.", format(n_snps, big.mark = ",")),
  "G-MoQ Stage 2 is ordinary least squares with the pooled Stage 1 covariance (Methods).",
  "**a**, Manhattan plots of the G-MoQ joint Hotelling test of the three mixture parameters (above the axis) and",
  "REGENIE on log MFI (below), for MCV VP1 (top) and JCV VP1 (bottom). Grey points, the largest −log₁₀ P in each",
  "1-Mb bin; coloured points, individual variants with joint P < 10⁻⁵, shown on both sides. Dashed lines,",
  if (any(fg$capped)) sprintf("P = 5 × 10⁻⁸; values above %d are drawn at %d as triangles.", Y_CAP, Y_CAP)
  else "P = 5 × 10⁻⁸.",
  "**b**, Quantile–quantile plots of the joint test, the G-MoQ mean-effect test and REGENIE; each point is one",
  "quantile of the observed P-value distribution; dashed line, identity. Genomic inflation factors (λ, median):",
  paste0(lambda_txt("mcv_vp1_log_all"), "; ", lambda_txt("jcv_vp1_log_all"), "."),
  sprintf("Counting loci within 500 kb, the joint and mean tests and REGENIE find the same %d MCV VP1 and %d JCV VP1",
          n_loci[["mcv_vp1_log_all"]], n_loci[["jcv_vp1_log_all"]]),
  "loci, with none unique to either method. P values are unadjusted for multiple comparisons.",
  "JCV, JC polyomavirus; MCV, Merkel cell polyomavirus; MFI, mean fluorescence intensity; VP1, viral capsid protein 1.")
save_figure(fig, "figS09_antibody_manhattan_qq", source = src, caption = cap, width = 183, height = 170,
            results_used = "04_antibodies/decomposition")
