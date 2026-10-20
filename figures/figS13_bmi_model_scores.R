#!/usr/bin/env Rscript
# Supplementary Figure S13: per-lead model selection for BMI lead variants
#
# What: each of the 639 BMI lead profiles (log scale) is fitted with seven Stage 2 models and
#   scored by an AIC-like criterion that uses the full Stage 1 covariance. For most leads the
#   runner-up model is within 2 units of the winner, so the winning model does not identify
#   the mechanism (compare Fig. 4).
#   a: table per winning model: wins, share, median gap to runner-up, share of gaps below 2,
#      and share of wins whose residual lack-of-fit test rejects.
#   b: distribution of the gap between winner and runner-up.
# Paper: Supplementary Fig. S13
# Input: results/05_bmi/model_scores/{lead_winners,lead_scores,model_definitions,gap_summary}.csv
# Output: paper_assets/figS13_bmi_model_scores.{png,pdf}, _source.csv, _caption.md
# Runs: locally from results/ in about 2 seconds

suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(data.table)})
source(file.path(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE), "figures", "theme.R"))

## 1. Load (primary lead set, log BMI, AIC-like criterion) -------------------------------------
winners <- as.data.table(read_result("05_bmi/model_scores", "lead_winners"))[set == "primary" & scale == "log" & criterion == "AIC"]
scores <- as.data.table(read_result("05_bmi/model_scores", "lead_scores"))[set == "primary" & scale == "log"]
models <- as.data.table(read_result("05_bmi/model_scores", "model_definitions"))[set == "primary" & scale == "log"]
gap_summary <- as.data.table(read_result("05_bmi/model_scores", "gap_summary"))[set == "primary" & scale == "log" & criterion == "AIC"]

## 2. Panel a: summary per winning model -------------------------------------------------
# Residual lack of fit of the winning model: Q_full against chi-square on 45 - k degrees of
# freedom (k, the model's effective number of parameters), one-sided.
won <- merge(winners[, .(snp_id, winner)], scores[model != "builder_shifted", .(snp_id, model, Q_full)], by = "snp_id")
won <- won[model == winner]
won <- merge(won, models[, .(model, k)], by = "model")
won[, lof_rejected := pchisq(Q_full, 45 - k, lower.tail = FALSE) < 0.05]
stopifnot(nrow(won) == nrow(winners))
tab <- won[, .(
  n = .N,
  pct = 100 * .N / nrow(winners),
  median_gap = median(winners[winner == unique(model), gap]),
  pct_gap_lt_2 = 100 * mean(winners[winner == unique(model), gap < 2]),
  pct_lof_rejected = 100 * mean(lof_rejected)
), by = model][order(-n)]
model_labels <- c(additive = "Additive", log_additive = "Log-additive", shifted = "Log-additive, shifted baseline",
                  both = "Additive + log-additive", logscale_vqtl = "vQTL (log scale)",
                  shifted_lognormal_full = "Additive + vQTL", three_components = "Shift + spread + bend")
tab[, model_name := model_labels[model]]
tab[, `:=`(n_label = as.character(n), pct_label = sprintf("%.1f%%", pct), gap_label = sprintf("%.2f", median_gap),
           gap2_label = sprintf("%.1f%%", pct_gap_lt_2), lof_label = sprintf("%.1f%%", pct_lof_rejected))]
# Top row: all leads together.
tab <- rbind(tab, data.table(
  model = "TOTAL", n = nrow(winners), pct = 100.0, median_gap = median(winners$gap),
  pct_gap_lt_2 = 100 * mean(winners$gap < 2), pct_lof_rejected = 100 * mean(won$lof_rejected),
  model_name = "All leads", n_label = as.character(nrow(winners)), pct_label = "100%",
  gap_label = sprintf("%.2f", median(winners$gap)), gap2_label = sprintf("%.1f%%", 100 * mean(winners$gap < 2)),
  lof_label = sprintf("%.1f%%", 100 * mean(won$lof_rejected))), fill = TRUE)
tab[, model_f := factor(model_name, levels = c("All leads", rev(tab$model_name[tab$model != "TOTAL"])))]
# Drawn as a table of tiles; short column headers are defined in the legend.
hdr <- c(n_label = "Wins", pct_label = "%", gap_label = "Median gap", gap2_label = "Gap < 2",
         lof_label = "Lack of fit")
tab_long <- melt(tab, id.vars = c("model_f", "model"), measure.vars = names(hdr),
                 variable.name = "field", value.name = "label")
tab_long[, field := factor(hdr[as.character(field)], levels = unname(hdr))]
pA <- ggplot(tab_long, aes(field, model_f)) +
  geom_tile(data = tab_long[model != "TOTAL"], fill = paper_palette[["shift"]], alpha = 0.12,
            colour = "white", linewidth = 0.6) +
  geom_tile(data = tab_long[model == "TOTAL"], fill = NA, colour = "grey50", linewidth = 0.4) +
  geom_text(data = tab_long[model != "TOTAL"], aes(label = label), size = 6 / .pt, colour = "grey20") +
  geom_text(data = tab_long[model == "TOTAL"], aes(label = label), size = 6 / .pt, colour = "grey20",
            fontface = "bold") +
  scale_x_discrete(position = "top", expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL, tag = "a") +
  theme_paper() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(colour = "grey20"), axis.text.y = element_text(colour = "grey20"))

## 3. Panel b: distribution of the gap -----------------------------------------------------
pB <- ggplot(winners, aes(gap)) +
  geom_histogram(binwidth = 0.1, boundary = 0, fill = paper_palette[["shift"]], colour = "white", linewidth = 0.1) +
  geom_vline(xintercept = 2, linetype = 2, colour = "grey30", linewidth = 0.35) +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.03))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "AIC-like gap to runner-up model", y = "Lead variants", tag = "b") +
  theme_paper()

## 4. Assemble and save ---------------------------------------------------------------------
fig <- (pA | pB) + plot_layout(widths = c(1.35, 1))
src <- rbind(
  tab[, .(panel = "A", item = model, n, pct, median_gap, pct_gap_lt_2, pct_lof_rejected)],
  winners[, .(panel = "B", item = "gap_dist", n = .N, pct = NA_real_, median_gap = median(gap), pct_gap_lt_2 = NA_real_, pct_lof_rejected = NA_real_)])
n_lead <- nrow(winners); n_within <- sum(winners$gap < 2)
stopifnot(isTRUE(all.equal(gap_summary$frac_gap_lt_2, n_within / n_lead)))
cap <- paste(
  "Figure S13. Per-lead model selection for BMI lead variants.",
  sprintf("Each of the %d BMI lead variants was fitted on the log scale with %s models of its Stage 1 quantile-effect",
          n_lead, c("one","two","three","four","five","six","seven","eight","nine")[nrow(models)]),
  "profile, by ordinary least squares with the pooled Stage 1 covariance, and scored by an AIC-like criterion using",
  "the full covariance (Methods).",
  "**a**, Per winning model: the number of leads it wins (Wins) and their share (%); the median AIC-like gap from the",
  "winner to the runner-up (Median gap); the percentage of its wins with a gap below 2 (Gap < 2); and the",
  "percentage of its wins whose residual lack-of-fit test rejects at P < 0.05 (Lack of fit; one-sided χ² test on",
  "45 − k degrees of freedom, k the model's effective number of parameters; not adjusted for multiple comparisons).",
  "Top row, all leads.",
  sprintf("**b**, Distribution of the gap across the %d leads; bars, leads per 0.1-unit bin; dashed line, a gap of 2.",
          n_lead),
  sprintf("%.1f%% of leads (%d of %d) have a runner-up within 2 units, so the winning model need not identify the",
          100 * n_within / n_lead, n_within, n_lead),
  "mechanism (Fig. 4b); near-ties mostly pair a model with its own one-parameter extension.",
  "AIC, Akaike information criterion; BMI, body mass index; vQTL, variance quantitative trait locus."
)
save_figure(fig, "figS13_bmi_model_scores", src, cap, width = 183, height = 110,
            results_used = "05_bmi/model_scores")
