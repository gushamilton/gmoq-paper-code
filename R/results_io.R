# Reading settings and aggregate results, and writing results safely.
#
# Every analysis script reads its settings from config.yaml and writes its outputs as CSV
# files under results/<section>/<analysis>/. Figures and tables read only these files.

# Root of this repository: $GMOQ_REPO if set, otherwise the git top level, otherwise the
# parent of R/.
repo_path <- function(...) {
  root <- Sys.getenv("GMOQ_REPO", "")
  if (!nzchar(root)) {
    root <- tryCatch(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE,
                             stderr = FALSE), error = function(e) "", warning = function(w) "")
    if (length(root) != 1L || !nzchar(root)) root <- normalizePath(file.path(.gmoq_r_dir, ".."))
  }
  file.path(normalizePath(root), ...)
}
repo_root <- function() repo_path()

# Settings: paths to the UK Biobank-derived inputs and the analysis choices.
# YAML would otherwise read the bare keys y / n as booleans, so they are kept as strings.
read_config <- function(path = Sys.getenv("GMOQ_CONFIG", repo_path("config.yaml"))) {
  keep_yn <- function(default) function(x) if (x %in% c("y", "Y", "n", "N")) x else default
  yaml::read_yaml(path, handlers = list("bool#yes" = keep_yn(TRUE), "bool#no" = keep_yn(FALSE)))
}

# The 45 quantile levels used throughout: 0.02, 0.0418, ..., 0.98.
tau_grid <- function(cfg = NULL) {
  g <- if (is.null(cfg)) list(tau_from = 0.02, tau_to = 0.98, n_tau = 45) else cfg$grid
  seq(g$tau_from, g$tau_to, length.out = g$n_tau)
}

# Which Stage 2 weight design is reported for a mixture fit: 'model' (every weight column
# divided by the fitted mixture density) everywhere except unshrunk testosterone, which uses
# 'hybrid' (the membership column divided by the Stage 1 kernel density instead).
mixture_design <- function(cfg, fit) {
  mw <- cfg$estimation$mixture_weights
  f <- gsub(".", "_", as.character(fit), fixed = TRUE)
  out <- rep(mw$default, length(f))
  overrides <- unlist(mw$overrides)
  hit <- f %in% names(overrides)
  out[hit] <- overrides[f[hit]]
  unname(out)
}

# Read results/<analysis>/<name>.csv, e.g. read_result("06_traits/moments", "lead_moments").
read_result <- function(analysis, name, root = repo_path("results")) {
  utils::read.csv(file.path(root, analysis, paste0(name, ".csv")),
                  stringsAsFactors = FALSE, check.names = FALSE)
}

# Write a data frame to <out>/<analysis>/<name>.csv. The output root is $GMOQ_OUT if set
# (useful for checking a rerun against the published results), otherwise results/.
# Applies a standard disclosure check before writing: no identifier-like columns and no count
# cells between 1 and 9.
write_result <- function(df, analysis, name, root = Sys.getenv("GMOQ_OUT", repo_path("results")),
                         count_cols = c("count", "n_cell", "cell_count")) {
  id_like <- grepl("^(eid|iid|fid|sample|participant|id)$|_eid$|_iid$", tolower(names(df)))
  if (any(id_like)) stop("refusing identifier-like column(s): ", paste(names(df)[id_like], collapse = ", "))
  for (cc in intersect(count_cols, names(df))) {
    v <- df[[cc]]
    if (any(!is.na(v) & v > 0 & v < 10)) stop("refusing counts below 10 in column '", cc, "'")
  }
  dir <- file.path(root, analysis)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(dir, paste0(name, ".csv"))
  utils::write.csv(df, path, row.names = FALSE, na = "")   # missing values as empty fields
  invisible(path)
}
