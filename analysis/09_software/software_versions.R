#!/usr/bin/env Rscript
# Software versions
#
# What: lists the software and versions used for the paper's results: the two G-MoQ
#   packages (Stage 1 in Python, Stage 2 in R, with the releases named in config.yaml),
#   the R and Python environments, and the external tools (PLINK, REGENIE, LDSC,
#   mixsmsn, sn). The versions are those of the cluster environment in which the analyses
#   ran, listed below. Run with --check on that cluster to compare the list against the
#   live environment.
# Paper: Supplementary Table S8
# Input: config.yaml (package releases and the Stage 1 Python environment)
# Output: results/09_software/versions/{package_pins,software}.csv
# Runs: locally in about 1 second (with --check: on the cluster where the analyses ran)

.repo <- Sys.getenv("GMOQ_REPO", system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE))
source(file.path(.repo, "R", "load.R"))
cfg <- read_config()
check_environment <- "--check" %in% commandArgs(trailingOnly = TRUE)

## 1. The G-MoQ packages
# Stage 1 produced every quantile GWAS. Stage 2 is the current release; the earlier release
# produced some analyses, and its default output is identical to the current one.
pk <- cfg$package
earlier <- pk$stage2_history[[1]]
package_pins <- data.table(
  role = c("stage1", "stage2 (current pin)", paste0("stage2 history (", earlier$short, ")")),
  commit = c(pk$stage1$commit, pk$stage2$commit, earlier$commit),
  short = c(pk$stage1$short, pk$stage2$short, earlier$short),
  version = c(NA, pk$stage2$r_package_version, earlier$r_package_version),
  lib = c(pk$stage1$root, pk$stage2$lib, earlier$lib),
  pin_json_commit = NA_character_, pin_json_matches = NA, install_job = NA_character_,
  note = c("Stage 1 quantile regressions (RIF-OLS) for every GWAS in the paper",
           "gmoq: Stage 2 projections, mixture models and lack-of-fit tests",
           "gmoq, earlier release: model-profile angles, BMI polygenic-score interactions and the weight-design simulation; default output identical to the current release"))

## 2. Versions in the analysis environment
software <- data.table(
  component = c("R", "data.table", "yaml", "jsonlite", "CompQuadForm", "testthat",
                "stage1 env: python", "stage1 env: numpy", "stage1 env: scipy",
                "stage1 env: pandas", "stage1 env: gmoq",
                "plink2 (lead clumping)", "Stage 1 density", "Stage 1 blocks / taus",
                "mixsmsn", "sn", "REGENIE", "LDSC"),
  version = c("R version 4.5.1 (2025-06-13)", "1.18.2.1", "2.3.12", "2.0.0", "1.4.4", "3.2.3",
              "3.10.19", "2.2.6", "1.15.2", "2.3.3", "0.1.0",
              "plink2 2.0.0-a.6.8LM", "Silverman Gaussian KDE (--density-method)",
              paste(cfg$grid$n_tau, "taus on", cfg$grid$tau_from, "-", cfg$grid$tau_to),
              "1.1.11", "2.1.1", "4.1.2", "2.0.0"),
  source = c(rep("R on the analysis cluster", 6), rep("Stage 1 Python environment", 5),
             "lead-variant clumping log", "Stage 1 default", "config.yaml grid",
             "R on the analysis cluster", "R on the analysis cluster", "REGENIE log", "LDSC log"))

## 3. Optional: compare the list with the live environment
if (check_environment) {
  live <- c(R = R.version.string,
            sapply(c("data.table", "yaml", "jsonlite", "CompQuadForm", "testthat", "mixsmsn", "sn"),
                   function(p) tryCatch(as.character(packageVersion(p)), error = function(e) NA)))
  py <- system2(pk$stage1$python, c("-c", shQuote(paste(
    "import sys, importlib; print('python|' + sys.version.split()[0])",
    "for m in ['numpy', 'scipy', 'pandas', 'gmoq']:",
    "    print('stage1 env: ' + m + '|' + str(getattr(importlib.import_module(m), '__version__', 'importable')))",
    sep = "\n"))), stdout = TRUE)
  py <- do.call(rbind, strsplit(py, "|", fixed = TRUE))
  live <- c(live, setNames(py[, 2], sub("^python$", "stage1 env: python", py[, 1])))
  listed <- setNames(software$version, software$component)[names(live)]
  differ <- names(live)[!is.na(live) & live != listed]
  if (length(differ)) warning("listed version differs from the environment for: ", paste(differ, collapse = ", "))
}

write_result(package_pins, "09_software/versions", "package_pins")
write_result(software, "09_software/versions", "software")
