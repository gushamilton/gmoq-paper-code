# Shared look and helpers for the paper's figures and supplementary tables
#
# What: every script in figures/ and tables/ sources this file. It loads the repository's
#   shared functions (R/load.R: read_result(), read_config(), mixture_design(), tau_grid()),
#   defines the journal figure style (Nature Portfolio artwork rules: 5-7 pt text, bold
#   lower-case panel letters in the margin, no gridlines, unboxed keys), a few plotting
#   helpers, and save_figure() / save_table(), which write each item and its legend to
#   paper_assets/.
# Paper: all main and supplementary figures and supplementary tables
# Input: nothing by itself; the scripts read aggregate CSV files from results/.
# Output: none by itself; scripts write paper_assets/<id>.png, <id>.pdf, <id>_source.csv
#   (figures) or <id>.csv (tables), and <id>_caption.md.
# Runs: locally; sourced at the top of every script in figures/ and tables/

suppressPackageStartupMessages(library(ggplot2))

repo <- local({
  out <- tryCatch(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE, stderr = FALSE),
                  error = function(e) character(0))
  if (length(out)) out else normalizePath(".")
})
source(file.path(repo, "R", "load.R"))

## 1. Colours and theme -------------------------------------------------------------------

# One colour per recurring idea, so the same thing has the same colour in every figure.
paper_palette <- c(
  shift = "#1b6ca8", spread = "#c0392b", bend = "#7d5ba6", residual = "#6b6b6b",
  female = "#b5446e", male = "#2a7f62", mixture = "#d98b2b", regenie = "#444444",
  gmoq = "#1b6ca8", null = "#9e9e9e", highlight = "#e4572e")

# Base theme: 7 pt text, axis lines and ticks, no gridlines, bare strip labels, no panel
# titles (the legend says what each panel shows), panel letters 8 pt bold in the margin.
theme_paper <- function(base_size = 7) {
  theme_minimal(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.line = element_line(linewidth = 0.3, colour = "grey25", lineend = "round"),
          axis.ticks = element_line(linewidth = 0.3, colour = "grey25"),
          axis.ticks.length = unit(1.5, "pt"),
          axis.text = element_text(colour = "grey20"),
          strip.text = element_text(hjust = 0.5, size = rel(1)),
          strip.background = element_blank(),
          plot.title = element_blank(),
          plot.subtitle = element_blank(),
          plot.tag = element_text(face = "bold", size = 8, hjust = 0, vjust = 1),
          plot.tag.position = "topleft", plot.tag.location = "margin",
          plot.margin = margin(9, 6, 4, 4),
          legend.position = "bottom", legend.key.size = unit(0.75, "lines"),
          legend.text = element_text(size = rel(0.9)), legend.title = element_text(size = rel(0.9)),
          legend.margin = margin(1, 2, 1, 2), legend.box.spacing = unit(2, "pt"))
}

# Put the key inside the plotting area: (x, y) is the anchor in panel coordinates (0-1) and
# `just` says which corner of the key sits on the anchor. The key has no box and no title.
legend_inside <- function(x = 0.98, y = 0.98, just = c(1, 1)) {
  theme(legend.position = "inside", legend.position.inside = c(x, y), legend.justification = just,
        legend.background = element_blank(),
        legend.key = element_blank(), legend.title = element_blank(),
        legend.key.height = unit(0.65, "lines"), legend.spacing.y = unit(0, "pt"))
}

## 2. Smoothing for display ---------------------------------------------------------------

# Smooth curve through tabulated points (for example the 45 quantile levels), for drawing
# only. Monotone Hermite interpolation (Fritsch-Carlson) passes through every point and
# never overshoots between neighbours, so it removes the kinks of a polyline without
# inventing wiggles. Returns data.table(x, y) on an n-point grid, separately for each
# combination of the `by` columns (which are carried along).
smooth_curve <- function(d, x, y, by = character(0), n = 400) {
  d <- data.table::as.data.table(d)
  interpolate <- function(xx, yy) {
    o <- order(xx); xx <- xx[o]; yy <- yy[o]
    ok <- is.finite(xx) & is.finite(yy); xx <- xx[ok]; yy <- yy[ok]
    if (length(xx) < 3) return(list(x = xx, y = yy))
    g <- seq(min(xx), max(xx), length.out = n)
    list(x = g, y = stats::splinefun(xx, yy, method = "monoH.FC")(g))
  }
  if (!length(by)) return(data.table::as.data.table(interpolate(d[[x]], d[[y]])))
  d[, interpolate(get(x), get(y)), by = by]
}

# Smooth density from a histogram with bins [lower, upper) and heights `density`: the
# histogram convolved with a Gaussian kernel of bandwidth bw, evaluated on `grid`. It keeps
# the area of the histogram, so smoothed parts still add up to the smoothed whole.
smooth_hist_density <- function(lower, upper, density, grid, bw) {
  ok <- is.finite(lower) & is.finite(upper) & is.finite(density)
  lower <- lower[ok]; upper <- upper[ok]; density <- density[ok]
  vapply(grid, function(x) sum(density * (stats::pnorm((upper - x) / bw) - stats::pnorm((lower - x) / bw))), 0)
}

## 3. Mixture weight design ---------------------------------------------------------------

# Results for the two-component mixtures are stored under both Stage 2 weight designs
# ('model': every weight column divided by the fitted mixture density; 'hybrid': the
# membership column divided by the Stage 1 kernel density instead). Keep only the rows of
# the design reported for each fit, as set in config.yaml. The fit is read from column
# `by`, or given once as `fit`. Tables without a weight_design column are returned as is.
pick_design <- function(x, by = "arm", fit = NULL, cfg = read_config()) {
  x <- data.table::as.data.table(x)
  if (!"weight_design" %in% names(x)) return(x)
  want <- if (!is.null(fit)) rep(mixture_design(cfg, fit), nrow(x)) else mixture_design(cfg, x[[by]])
  x[x$weight_design == want]
}

## 4. Single-variant layout (Figs 2d,e and 3b-d) ------------------------------------------

# Left half: forest plot of per-allele estimates with 95% confidence intervals.
#   fd: data.frame(label, est, se, ck, block). Rows are drawn top to bottom in the order
#       given; `ck` picks the colour from `cols`; alternate blocks are shaded.
snp_forest <- function(fd, cols, xlab = "Per-allele effect (95% CI)", block_fill = "grey95") {
  fd <- data.table::as.data.table(fd)
  fd[, y := rev(seq_len(.N))]
  fd[, `:=`(lo = est - 1.96 * se, hi = est + 1.96 * se)]
  bands <- fd[, .(ymin = min(y) - 0.5, ymax = max(y) + 0.5), by = block]
  bands[, shade := seq_len(.N) %% 2 == 1]
  ggplot(fd) +
    geom_rect(data = bands[shade == TRUE], aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
              fill = block_fill, inherit.aes = FALSE) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey55") +
    geom_errorbar(aes(y = y, xmin = lo, xmax = hi, colour = ck), width = 0, linewidth = 0.55, orientation = "y") +
    geom_point(aes(est, y, colour = ck), size = 1.6) +
    scale_y_continuous(breaks = fd$y, labels = fd$label, expand = expansion(add = 0.5)) +
    scale_colour_manual(values = cols, guide = "none") +
    labs(x = xlab, y = NULL) +
    theme_paper() + theme(panel.grid.major.y = element_blank(), axis.line.y = element_blank(),
                          axis.ticks.y = element_blank())
}

# Right half: the variant's Stage 1 quantile effects (grey points with 95% CI) and the
# Stage 2 fit as smooth curves: the total (series "total") and each parameter's
# contribution W[, p] * theta_p.
#   obs:  data.frame(tau, beta, se);  fits: data.frame(tau, series, value).
#   ylim: display range; points outside it are drawn at the edge as triangles.
#   headroom: extra space above the data (fraction of its range) for a key in the panel.
snp_profile <- function(obs, fits, cols, labels = NULL, ylim = NULL,
                        ylab = "Stage 1 quantile effect (95% CI)", xlab = "Quantile level tau",
                        legend_pos = c(0.98, 0.02), legend_just = c(1, 0), headroom = 0) {
  obs <- data.table::as.data.table(obs); fits <- data.table::as.data.table(fits)
  if (is.null(ylim) && headroom > 0) {
    r <- range(c(obs$beta - 1.96 * obs$se, obs$beta + 1.96 * obs$se, fits$value, 0), na.rm = TRUE)
    ylim <- c(r[1] - 0.03 * diff(r), r[2] + headroom * diff(r))
  }
  obs[, `:=`(lo = beta - 1.96 * se, hi = beta + 1.96 * se, clipped = FALSE)]
  if (!is.null(ylim)) {
    obs[, clipped := beta < ylim[1] | beta > ylim[2]]
    obs[, `:=`(beta = pmin(pmax(beta, ylim[1]), ylim[2]), lo = pmin(pmax(lo, ylim[1]), ylim[2]),
               hi = pmin(pmax(hi, ylim[1]), ylim[2]))]
  }
  sm <- smooth_curve(fits, "tau", "value", by = "series")
  sm[, series := factor(series, levels = names(cols)[names(cols) %in% unique(series)])]
  lty <- setNames(rep("solid", nlevels(sm$series)), levels(sm$series))
  lw <- setNames(ifelse(levels(sm$series) == "total", 0.8, 0.55), levels(sm$series))
  key_labels <- if (is.null(labels)) waiver() else labels
  ggplot() +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey55") +
    geom_linerange(data = obs, aes(tau, ymin = lo, ymax = hi), colour = "grey70", linewidth = 0.3) +
    geom_point(data = obs, aes(tau, beta, shape = clipped), colour = "grey35", size = 0.8) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
    geom_line(data = sm[series != "total"], aes(x, y, colour = series, linetype = series, linewidth = series),
              lineend = "round") +
    geom_line(data = sm[series == "total"], aes(x, y, colour = series, linetype = series, linewidth = series),
              lineend = "round") +
    scale_colour_manual(values = cols, labels = key_labels, name = NULL, breaks = levels(sm$series)) +
    scale_linetype_manual(values = lty, labels = key_labels, name = NULL, breaks = levels(sm$series)) +
    scale_linewidth_manual(values = lw, labels = key_labels, name = NULL, breaks = levels(sm$series)) +
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", ".25", ".5", ".75", "1")) +
    coord_cartesian(xlim = c(0, 1), ylim = ylim) +
    labs(x = xlab, y = ylab) +
    theme_paper() + legend_inside(legend_pos[1], legend_pos[2], legend_just)
}

## 5. Saving ------------------------------------------------------------------------------

assets_dir <- file.path(repo, "paper_assets")

# Refuse anything that looks like a participant identifier in a table leaving results/.
check_no_ids <- function(df) {
  id_like <- grepl("^(eid|iid|fid|sample|participant|id)$|_eid$|_iid$", tolower(names(df)))
  if (any(id_like)) stop("refusing identifier-like column(s): ", paste(names(df)[id_like], collapse = ", "))
}

# The legend file: a heading with the item id, the legend text, and one line naming the
# results folders the item was drawn from.
write_caption <- function(id, caption, results_used, dir) {
  src <- if (length(results_used)) paste0("Source data: ", paste0("results/", results_used, collapse = ", "), ".")
  writeLines(c(paste0("# ", id), "", caption, if (length(src)) c("", src)),
             file.path(dir, paste0(id, "_caption.md")))
}

# Write paper_assets/<id>.png (450 dpi), <id>.pdf, <id>_source.csv (the numbers plotted)
# and <id>_caption.md. Width and height in mm (journal widths: 89, 120, 136 or 183 mm).
#   results_used: the results folders read, e.g. "03_testosterone/recovery".
save_figure <- function(plot, id, source, caption, width = 183, height = 120,
                        results_used = character(0), dir = assets_dir, dpi = 450) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  check_no_ids(source)
  png <- file.path(dir, paste0(id, ".png")); pdf <- file.path(dir, paste0(id, ".pdf"))
  ggsave(png, plot, width = width, height = height, units = "mm", dpi = dpi, bg = "white",
         device = ragg::agg_png)
  # cairo_pdf needs X11 libraries that some Macs lack; fall back to the base PDF device.
  tryCatch(ggsave(pdf, plot, width = width, height = height, units = "mm", device = grDevices::cairo_pdf),
           error = function(e)
             ggsave(pdf, plot, width = width, height = height, units = "mm", device = grDevices::pdf),
           warning = function(w)
             ggsave(pdf, plot, width = width, height = height, units = "mm", device = grDevices::pdf))
  utils::write.csv(source, file.path(dir, paste0(id, "_source.csv")), row.names = FALSE, na = "")
  write_caption(id, caption, results_used, dir)
  invisible(c(png = png, pdf = pdf))
}

# Write paper_assets/<id>.csv and <id>_caption.md for a supplementary table.
save_table <- function(df, id, caption, results_used = character(0), dir = assets_dir) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  check_no_ids(df)
  utils::write.csv(df, file.path(dir, paste0(id, ".csv")), row.names = FALSE, na = "")
  write_caption(id, caption, results_used, dir)
  invisible(file.path(dir, paste0(id, ".csv")))
}
