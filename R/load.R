# Load the packages and shared functions used by every script in this repository:
#   source(file.path(<repo>, "R", "load.R"))
suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
})
.gmoq_r_dir <- local({
  here <- NULL
  for (i in rev(seq_len(sys.nframe()))) {
    of <- sys.frame(i)$ofile
    if (!is.null(of) && basename(of) == "load.R") { here <- dirname(normalizePath(of)); break }
  }
  if (is.null(here)) here <- file.path(Sys.getenv("GMOQ_REPO", getwd()), "R")
  here
})
# Every other R/*.R file: results I/O and the statistics helpers, one topic per file.
for (.f in setdiff(list.files(.gmoq_r_dir, pattern = "\\.R$"), "load.R")) {
  source(file.path(.gmoq_r_dir, .f), local = FALSE)
}
rm(.f)
