#!/usr/bin/env Rscript
# Supplementary Table S8: software and versions
#
# What: one row per software package used for the results in the paper, with its version and
#   what it was used for. Versions are read from what was recorded in the analysis
#   environment.
# Paper: Supplementary Table S8
# Input: results/09_software/versions/{software,package_pins}.csv
# Output: paper_assets/tableS8_software.csv, _caption.md
# Runs: locally from results/ in about 1 second

suppressPackageStartupMessages({library(data.table)})
source(file.path(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE), "figures", "theme.R"))

## 1. Versions recorded in the software results ----------------------------------------------
software <- as.data.table(read_result("09_software/versions", "software"))
pins <- as.data.table(read_result("09_software/versions", "package_pins"))
version_of <- function(component_name) {
  v <- software[component == component_name, version]; stopifnot(length(v) == 1, nzchar(v)); v
}
stage2_current <- pins[role == "stage2 (current pin)"]
stage2_earlier <- pins[startsWith(role, "stage2 history")]
stopifnot(nrow(stage2_current) == 1, nrow(stage2_earlier) == 1)
r_version <- sub("^R version ([0-9.]+).*", "\\1", version_of("R"))
plink_version <- sub("^plink2 ", "", version_of("plink2 (lead clumping)"))

## 2. Tool versions recorded from the analysis logs ------------------------------------------
regenie_version <- version_of("REGENIE"); ldsc_version <- version_of("LDSC")

## 3. Table ---------------------------------------------------------------------------------
out <- data.table(
  component = c("G-MoQ Stage 1 (Python package)", "G-MoQ Stage 2 (R package)",
                "G-MoQ Stage 2 (R package), earlier release",
                "Python", "NumPy", "SciPy", "pandas",
                "R", "data.table", "CompQuadForm", "mixsmsn", "sn",
                "PLINK", "REGENIE", "LDSC"),
  version = c(version_of("stage1 env: gmoq"), stage2_current$version, stage2_earlier$version,
              version_of("stage1 env: python"), version_of("stage1 env: numpy"), version_of("stage1 env: scipy"),
              version_of("stage1 env: pandas"),
              r_version, version_of("data.table"), version_of("CompQuadForm"), version_of("mixsmsn"), version_of("sn"),
              plink_version, regenie_version, ldsc_version),
  use = c("Stage 1 quantile regressions (all traits)",
          "Stage 2 projections, mixture models, lack of fit",
          "Model-profile angles (Fig. 4b), BMI polygenic-score interactions (Fig. 6, Supplementary Fig. S14) and the weight-design simulation (Supplementary Fig. S7); default output identical to the current release",
          "Stage 1 environment", "Stage 1 environment", "Stage 1 environment", "Stage 1 environment",
          "Stage 2 and all downstream analyses", "Data handling", "Generalised chi-square (Davies) P values",
          "Skew-normal mixture fits (mixture densities, simulations)", "Skew-normal densities (simulations)",
          "Lead-variant clumping, LD, genotype extraction and the BMI polygenic score (clumping and thresholding)",
          "Reference GWAS", "Heritability and genetic correlation"))
stopifnot(!anyNA(out$version), all(nzchar(out$version)))

## 4. Save ----------------------------------------------------------------------------------
cap <- paste(
  "Table S8. Software and versions. One row per software package used for the results reported in the paper, with",
  "its release and the analyses it served (Methods). The earlier G-MoQ Stage 2 release produced the analyses listed",
  "against it; its default output is identical to the current release. Versions are those recorded in the analysis",
  "environment. Test, file-format and",
  "plotting software are not listed. BMI, body mass index; GWAS, genome-wide association study; LD, linkage",
  "disequilibrium; LDSC, LD score regression.")
save_table(as.data.frame(out), "tableS8_software", cap,
           results_used = "09_software/versions")
