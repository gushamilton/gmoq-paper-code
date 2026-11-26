#!/usr/bin/env Rscript
# Fig 1: G-MoQ conceptual overview (a-c) + reading a quantile profile through a mixture model (d-e).
#   (a) Stage 1: additive per-allele quantile shifts across AA (0), AG (1), GG (2);
#   (b) Stage 1 output: SNP x quantile summary-statistics table beta_hat(tau) (run once);
#   (c) Stage 2: model design matrix W (45 x p) converting beta_hat(tau) into parameters theta_hat;
#   (d) UKB testosterone distributions by sex across the shrinkage series (n = 399,710);
#   (e) Stage 2 weights W[tau, p] for female location, male location, and mixing parameter gamma.
suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(data.table)})
.repo <- system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE)
source(file.path(.repo, "figures", "theme.R"))

# ---- TOP CONCEPTUAL STRIP (a, b, c) ----------------------------------------------------------
col_aa <- "#8e9aab"; col_ag <- "#5899c4"; col_gg <- paper_palette[["gmoq"]]
theme_c <- function(base_size = 6.8) {
  theme_classic(base_size = base_size, base_family = "Helvetica") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "grey25"),
      axis.ticks = element_line(linewidth = 0.3, colour = "grey25"),
      axis.text = element_text(colour = "grey20", size = 5.6),
      axis.title = element_text(colour = "grey10", size = 6.4),
      legend.background = element_blank(), legend.key = element_blank(),
      legend.text = element_text(size = 5.5), legend.title = element_blank(),
      plot.tag = element_text(face = "bold", size = 8, colour = "grey10"),
      plot.margin = margin(3, 5, 3, 4)
    )
}

# a: AA / AG / GG additive per-allele quantile shifts
xg <- seq(-3.0, 4.4, length.out = 400)
da <- rbind(
  data.table(x = xg, y = dnorm(xg, 0.00, 0.95), geno = "AA (0)"),
  data.table(x = xg, y = dnorm(xg, 0.45, 1.09), geno = "AG (1)"),
  data.table(x = xg, y = dnorm(xg, 0.90, 1.23), geno = "GG (2)")
)
da[, geno := factor(geno, levels = c("AA (0)", "AG (1)", "GG (2)"))]
q_taus <- c(0.15, 0.50, 0.85)
q_aa <- qnorm(q_taus, 0.00, 0.95); q_ag <- qnorm(q_taus, 0.45, 1.09); q_gg <- qnorm(q_taus, 0.90, 1.23)
q_dt <- data.table(x0 = q_aa, x1 = q_ag, x2 = q_gg, y = c(0.14, 0.31, 0.14),
                   lbl = c("tau == 0.15", "tau == 0.50", "tau == 0.85"))
p0_a <- ggplot(da, aes(x, y, colour = geno, fill = geno)) +
  geom_area(alpha = 0.14, position = "identity", linewidth = 0.5) +
  geom_vline(xintercept = q_aa, linetype = "dotted", linewidth = 0.28, colour = col_aa) +
  geom_segment(data = q_dt, aes(x = x0, xend = x1, y = y, yend = y), inherit.aes = FALSE,
               colour = "grey20", linewidth = 0.42, arrow = arrow(length = unit(1.1, "mm"), type = "closed")) +
  geom_segment(data = q_dt, aes(x = x1, xend = x2, y = y, yend = y), inherit.aes = FALSE,
               colour = "grey20", linewidth = 0.42, arrow = arrow(length = unit(1.1, "mm"), type = "closed")) +
  geom_text(data = q_dt, aes(x = x1, y = y + 0.038, label = lbl), inherit.aes = FALSE,
            parse = TRUE, size = 1.8, colour = "grey15") +
  scale_colour_manual(values = c("AA (0)" = col_aa, "AG (1)" = col_ag, "GG (2)" = col_gg)) +
  scale_fill_manual(values = c("AA (0)" = col_aa, "AG (1)" = col_ag, "GG (2)" = col_gg)) +
  scale_y_continuous(limits = c(0, 0.48), expand = expansion(mult = c(0, 0))) +
  labs(x = "Phenotype Y", y = "Density", tag = "a") +
  theme_c() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                    legend.position = c(0.79, 0.85), legend.key.size = unit(2.3, "mm"))

# b: Stage 1 output table (SNPs x 45 quantiles)
snps <- c("SNP 1\n(shift)", "SNP 2\n(spread)", "SNP 3\n(mixture)", "SNP 4\n(LoD)", "...")
q_cols <- c("0.05", "0.25", "0.50", "0.75", "0.95")
vals <- matrix(c(
   0.04,  0.04,  0.04,  0.04,  0.04,
  -0.05, -0.02,  0.00,  0.02,  0.05,
   0.00,  0.00,  0.01,  0.05,  0.06,
   0.00,  0.00,  0.04,  0.04,  0.04,
     NA,    NA,    NA,    NA,    NA
), nrow = 5, byrow = TRUE)
tbl_dt <- CJ(r = 1:5, c = 1:5)
tbl_dt[, `:=`(snp = factor(snps[r], levels = rev(snps)),
              q   = factor(q_cols[c], levels = q_cols),
              val = mapply(function(i, j) vals[i, j], r, c))]
tbl_dt[, lbl := ifelse(is.na(val), "...", sprintf("%+.2f", val))]
p0_b <- ggplot(tbl_dt, aes(q, snp, fill = val)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = lbl, colour = !is.na(val) & abs(val) > 0.042), size = 1.85) +
  scale_colour_manual(values = c(`FALSE` = "grey15", `TRUE` = "white"), guide = "none") +
  scale_fill_gradient2(low = "#c4513b", mid = "#f5f6f7", high = paper_palette[["gmoq"]],
                       midpoint = 0, limits = c(-0.06, 0.06), na.value = "grey95", guide = "none") +
  scale_x_discrete(position = "top") +
  labs(x = expression("Stage 1 output: " * hat(beta)(tau) * " across 45 quantiles (run once)"),
       y = NULL, tag = "b") +
  theme_c() + theme(axis.line = element_blank(), axis.ticks = element_blank(),
                    axis.text.x.top = element_text(size = 5.6, face = "bold"),
                    axis.text.y = element_text(size = 5.5, lineheight = 0.9),
                    axis.title.x.top = element_text(size = 6.2, margin = margin(b = 3)))

# c: Stage 2 design matrix W (45 x p)
w_cols <- c("Female\nloc.", "Male\nloc.", "Mixing\n(gamma)")
w_rows <- c("0.05", "0.25", "0.50", "0.75", "0.95")
w_mat <- matrix(c(1.00, 0.00, -0.05, 1.00, 0.00, -0.15, 0.50, 0.50, -0.95, 0.00, 1.00, -0.15, 0.00, 1.00, -0.05),
                nrow = 5, byrow = TRUE)
w_dt <- CJ(r = 1:5, c = 1:3)
w_dt[, `:=`(q = factor(w_rows[r], levels = rev(w_rows)),
            par = factor(w_cols[c], levels = w_cols),
            val = mapply(function(i, j) w_mat[i, j], r, c))]
w_dt[, lbl := sprintf("%.2f", val)]
p0_c <- ggplot(w_dt, aes(par, q, fill = val)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = lbl, colour = abs(val) > 0.65), size = 1.85) +
  scale_colour_manual(values = c(`FALSE` = "grey15", `TRUE` = "white"), guide = "none") +
  scale_fill_gradient2(low = paper_palette[["mixture"]], mid = "#f5f6f7", high = paper_palette[["male"]],
                       midpoint = 0, limits = c(-1, 1), guide = "none") +
  scale_x_discrete(position = "top") +
  labs(x = expression("Stage 2 design matrix " * italic(W) * " (" * hat(beta) %~~% italic(W)*theta * ")"),
       y = expression("Quantile " * tau), tag = "c") +
  theme_c() + theme(axis.line = element_blank(), axis.ticks = element_blank(),
                    axis.text.x.top = element_text(size = 5.5, face = "bold", lineheight = 0.9),
                    axis.text.y = element_text(size = 5.5),
                    axis.title.x.top = element_text(size = 6.2, margin = margin(b = 3)))

pTop <- (p0_a | p0_b | p0_c) + plot_layout(widths = c(1.15, 1.15, 0.85))

# ---- TESTOSTERONE DENSITIES (d) AND STAGE 2 WEIGHTS (e) --------------------------------------
arm_lab <- c(testosterone = "No shrinkage", testosterone_shrink_0_2 = "Shrink 0.2",
             testosterone_shrink_0_5 = "Shrink 0.5", testosterone_shrink_0_7 = "Shrink 0.7")
h <- as.data.table(read_result("03_testosterone/recovery", "phenotype_histograms"))
h[, arm := gsub(".", "_", phenotype, fixed = TRUE)]
h <- h[is.finite(bin_lower) & is.finite(bin_upper) & !is.na(density)]
ng <- h[, .(n_group = n_group[1]), by = .(arm, group)]
nall <- ng[group == "all", setNames(n_group, arm)]
h[, share := n_group / nall[arm]]
h[, dens_plot := density * ifelse(group == "all", 1, share)]
h[, arm_f := factor(arm_lab[arm], levels = arm_lab)]
hs <- h[group != "all"]; ha <- h[group == "all"]
grid <- seq(-1.1, 3.4, length.out = 500)
sm <- h[, { bw <- 1.0 * median(bin_upper - bin_lower)
            .(x = grid, y = smooth_hist_density(bin_lower, bin_upper, dens_plot, grid, bw)) }, by = .(arm, arm_f, group)]
sm <- sm[y > 2e-3]
pD <- ggplot() +
  geom_area(data = sm[group != "all"], aes(x, y, fill = group), alpha = 0.5, position = "identity",
            colour = NA, outline.type = "full") +
  geom_line(data = sm[group != "all"], aes(x, y, colour = group), linewidth = 0.35) +
  geom_line(data = sm[group == "all"], aes(x, y, colour = group), linewidth = 0.6, lineend = "round") +
  facet_wrap(~arm_f, nrow = 1) +
  scale_fill_manual(values = c(female = paper_palette[["female"]], male = paper_palette[["male"]]),
                    labels = c(female = "Female", male = "Male"), name = NULL) +
  scale_colour_manual(values = c(all = "grey10", female = paper_palette[["female"]], male = paper_palette[["male"]]),
                      breaks = "all", labels = c(all = "Whole cohort"), name = NULL) +
  scale_x_continuous(breaks = -1:3) +
  guides(fill = guide_legend(order = 1), colour = guide_legend(order = 2)) +
  labs(x = "Testosterone (analysed scale)", y = "Density", tag = "d") +
  theme_paper() + legend_inside(0.02, 0.98, c(0, 1)) + theme(legend.key.width = unit(0.9, "lines"))

cfg <- read_config()
wall <- as.data.table(read_result("02_mixtures/stage2_mixtures", "weight_matrices"))[fit %in% names(arm_lab)]
w <- pick_design(wall, "fit", cfg = cfg)
D0 <- mixture_design(cfg, "testosterone"); Dalt <- setdiff(unlist(cfg$estimation$mixture_weight_designs), D0)[1]
arm_lab_d <- arm_lab
w[, arm_f := factor(arm_lab_d[fit], levels = arm_lab_d)]
wa <- wall[fit == "testosterone" & weight_design == Dalt & param == "gamma"]
wa[, arm_f := factor(arm_lab_d[fit], levels = arm_lab_d)]
w[, param_f := factor(param, levels = c("beta_left", "beta_right", "gamma"),
                      labels = c("Female location\n(beta_left)", "Male location\n(beta_right)", "Sex mixture\n(gamma)"))]
stopifnot(max(abs(w$W_max_chr - w$W_min_chr)) < 1e-6)
wm <- copy(w); wm[fit == "testosterone" & param == "gamma", keep := FALSE]
wm <- wm[is.na(keep)][, keep := NULL]
wg0 <- rbind(wa[, .(arm_f, param, tau, W_chr01, line = "model")],
             w[fit == "testosterone" & param == "gamma", .(arm_f, param, tau, W_chr01, line = "hybrid")])
curve_df <- function(d, extra = character(0)) smooth_curve(d, "tau", "W_chr01", by = c("arm_f", "param", extra), n = 600)
sl <- curve_df(wm[param != "gamma"])
sg <- curve_df(wm[param == "gamma"]); sg[, line := "model"]
sg0 <- curve_df(wg0, "line")
sg <- rbind(sg, sg0)
row_theme <- theme(strip.text.y.right = element_text(angle = 0, hjust = 0, face = "bold"))
xsc <- scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", ".25", ".5", ".75", "1"))
sl[, param_f := factor(param, levels = c("beta_left", "beta_right"), labels = c("Female\nlocation", "Male\nlocation"))]
pE1 <- ggplot(sl, aes(x, y, colour = param)) +
  geom_hline(yintercept = c(0, 1), linewidth = 0.25, linetype = 2, colour = "grey60") +
  geom_line(linewidth = 0.7, lineend = "round") +
  facet_grid(param_f ~ arm_f) + xsc +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  scale_colour_manual(values = c(beta_left = paper_palette[["female"]], beta_right = paper_palette[["male"]]), guide = "none") +
  labs(x = NULL, y = "Stage 2 weight", tag = "e") +
  theme_paper() + row_theme + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
sg[, param_f := "Mixing\n(gamma)"]
pE2 <- ggplot(sg, aes(x, y)) +
  geom_hline(yintercept = 0, linewidth = 0.25, linetype = 2, colour = "grey60") +
  geom_line(data = sg[line == "model"], colour = paper_palette[["mixture"]], linewidth = 0.7, lineend = "round") +
  geom_line(data = sg[line == "hybrid"], colour = "grey10", linewidth = 0.6, linetype = "22") +
  facet_grid(param_f ~ arm_f) + xsc +
  labs(x = "Quantile level tau", y = NULL) +
  theme_paper() + row_theme + theme(strip.text.x = element_blank())

pE <- pE1 / pE2 + plot_layout(heights = c(2, 1.15))
p <- (pTop / pD / pE) + plot_layout(heights = c(0.92, 0.95, 1.75))

src <- rbind(
  hs[, .(panel = "D", arm, series = group, x_lower = bin_lower, x_upper = bin_upper, value = dens_plot, n = n_group)],
  ha[, .(panel = "D", arm, series = "all", x_lower = bin_lower, x_upper = bin_upper, value = dens_plot, n = n_group)],
  w[, .(panel = "E", arm = fit, series = paste(param, weight_design), x_lower = tau, x_upper = tau, value = W_chr01, n = NA_real_)],
  wa[, .(panel = "E", arm = fit, series = paste(param, weight_design), x_lower = tau, x_upper = tau, value = W_chr01, n = NA_real_)])
cap <- paste(
  "Figure 1. Overview of the Genomic Method of Quantiles (G-MoQ) and application to a two-component testosterone mixture.",
  "**a**, Stage 1 estimates additive per-allele shifts across genotypes AA (0), AG (1) and GG (2) at each quantile tau.",
  "**b**, Stage 1 produces one summary-statistics table of SNP-by-quantile effects beta_hat(tau) across 45 quantiles.",
  "**c**, Stage 2 regresses each SNP's quantile profile on a model design matrix W (here female location, male location and mixing parameter gamma) to estimate mechanistic parameters theta.",
  "**d**, UK Biobank testosterone distributions (n = 399,710; 199,572 female, 200,138 male) in the unshrunk trait and three arms with male-female separation shrunk by 0.2, 0.5 and 0.7.",
  "**e**, Exact Stage 2 weights W[tau, p] for the fitted two-component skew-normal mixture across the four arms (solid, model-density weights; dashed in the unshrunk arm, hybrid membership weight; Methods).")

save_figure(p, "fig1_concept_testosterone", src, cap, width = 183, height = 205,
            results_used = c("03_testosterone/recovery", "02_mixtures/stage2_mixtures"))
