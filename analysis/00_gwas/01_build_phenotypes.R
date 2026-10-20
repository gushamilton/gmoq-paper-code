#!/usr/bin/env Rscript
# Build the phenotype, covariate and sample files for every GWAS in the paper
#
# What: reads UK Biobank fields, links them to the genetic identifiers and writes the
#   phenotype, covariate and keep (sample) files used by REGENIE and by the Stage 1
#   quantile GWAS: BMI and log BMI; height-adjusted weight with its 16 covariates;
#   log testosterone and its three shrunk arms; JC and Merkel cell polyomavirus VP1
#   antibody levels and serostatus; and the trait-survey phenotypes, including HbA1c
#   in participants without diabetes. Section 7 builds the samples that can only be
#   defined once the REGENIE step 1 predictions exist (run it after 02_regenie.sh step 1).
# Paper: Methods (phenotypes, samples and covariates); inputs to every figure and table
# Input: UK Biobank phenotype extract (tab-separated, columns f.<field>.<instance>.<array>),
#   the application-to-genetic-ID linker, the multiplex serology export, the imputed-data
#   sample file and, for section 7, the REGENIE step 1 LOCO files
# Output: participant-level phenotype, covariate and keep files in out_dir (stay on the
#   cluster; never copied into this repository)
# Runs: on a cluster with UK Biobank data
#
#   Rscript 01_build_phenotypes.R                 # sections 1-6
#   Rscript 01_build_phenotypes.R --after-step1   # section 7, after REGENIE step 1

suppressPackageStartupMessages(library(data.table))

# ---------------------------------------------------------------- user settings --
ukb_extract    <- "/path/to/ukb_phenotypes.tab.gz"      # UK Biobank main extract, f.<field>.0.0 columns
ukb_linker     <- "/path/to/linker_app_to_genetic.csv.gz" # columns "", ieu, app, eid (application eid -> genetic ID)
serology_file  <- "/path/to/antibody_data.tsv.gz"       # serology export, "<antigen> | Instance 0/1" columns
serology_linker <- "/path/to/linker_serology_app.csv"   # columns ieu, app for the serology application
bgen_sample    <- "/path/to/ukb_imputed/data.chr1-22.sample"
out_dir        <- "/path/to/gwas_inputs"
after_step1    <- "--after-step1" %in% commandArgs(trailingOnly = TRUE)

# REGENIE step 1 LOCO files (section 7 only); written by 02_regenie.sh
loco_log_bmi      <- file.path(out_dir, "regenie_step1/log_bmi/log_bmi_step1_1.loco")
loco_testosterone <- file.path(out_dir, "regenie_step1/testosterone/testosterone_step1_1.loco")
antibody_pred_list <- file.path(out_dir, "regenie_step1/antibody_qt_all/antibody_qt_all_step1_pred.list")

for (d in c("keep", "keep_panel")) dir.create(file.path(out_dir, d), recursive = TRUE, showWarnings = FALSE)
PCS <- paste0("PC", 1:10)
out <- function(...) file.path(out_dir, ...)
write_table <- function(x, path) fwrite(x, path, sep = " ", na = "NA", quote = FALSE)
write_keep <- function(fid, iid, path) fwrite(data.table(fid, iid), path, sep = " ", col.names = FALSE, quote = FALSE)

if (!after_step1) {

## 1. Core UK Biobank fields, linked to the genetic identifiers ---------------------
# Baseline (instance 0) values. Field IDs: 21003 age, 31 sex (0 female, 1 male), 50 standing
# height, 21002 weight, 21001 BMI, 30850 testosterone, 30710 CRP, 30870 triglycerides,
# 22009 genetic principal components 1-10.
core_fields <- c(age = "f.21003.0.0", sex = "f.31.0.0", height = "f.50.0.0", weight = "f.21002.0.0",
                 bmi = "f.21001.0.0", testosterone = "f.30850.0.0", crp = "f.30710.0.0",
                 triglycerides = "f.30870.0.0", setNames(paste0("f.22009.0.", 1:10), PCS))
raw <- fread(ukb_extract, select = c("f.eid", core_fields), sep = "\t")
setnames(raw, c("f.eid", core_fields), c("eid", names(core_fields)))
# Inner join on the application eid; the genetic identifier (ieu) becomes FID = IID.
linker <- unique(fread(ukb_linker, select = c("eid", "ieu")), by = "eid")
base <- merge(raw, linker, by = "eid")[!is.na(ieu)]
base[, `:=`(FID = ieu, IID = ieu)]
setcolorder(base, c("FID", "IID", "eid"))
base[, ieu := NULL]

## 2. BMI, weight and testosterone phenotypes -----------------------------------------
pheno <- copy(base)
pheno[, bmi := fifelse(is.finite(bmi) & bmi > 0, bmi, NA_real_)]
pheno[, log_bmi := fifelse(!is.na(bmi), log(bmi), NA_real_)]
pheno[, height := fifelse(is.finite(height) & height > 0, height, NA_real_)]
pheno[, weight := fifelse(is.finite(weight) & weight > 0, weight, NA_real_)]
# Testosterone is analysed on the log scale; values <= 0 are dropped.
pheno[, testosterone := fifelse(is.finite(testosterone) & testosterone > 0, log(testosterone), NA_real_)]
# Shrunk arms: male values are moved towards the female mean by a fraction F of the
# male-female difference in mean log testosterone, so the two sex components overlap more
# (F = 0.2, 0.5, 0.7). Female values are unchanged. Means use every participant with a value.
sex_shift <- pheno[sex == 1 & !is.na(testosterone), mean(testosterone)] -
             pheno[sex == 0 & !is.na(testosterone), mean(testosterone)]
for (shrink in c(0.2, 0.5, 0.7)) {
  arm <- paste0("testosterone_shrink_", shrink)
  pheno[, (arm) := testosterone]
  pheno[sex == 1 & !is.na(testosterone), (arm) := testosterone - shrink * sex_shift]
}
pheno_cols <- c("FID", "IID", "eid", "age", "sex", "height", PCS, "bmi", "log_bmi", "weight",
                "testosterone", paste0("testosterone_shrink_", c(0.2, 0.5, 0.7)))
pheno <- pheno[, ..pheno_cols]
write_table(pheno, out("phenotypes_paper.txt"))

## 3. Covariate files -------------------------------------------------------------------
# With sex: age, sex, PC1-PC10 (BMI, log BMI, trait survey, antibodies).
# Without sex: age, PC1-PC10 (testosterone, whose sex difference is the signal of interest,
# and the GWAS of sex itself).
cov_with_sex <- unique(pheno[, c("FID", "IID", "age", "sex", PCS), with = FALSE], by = c("FID", "IID"))
cov_no_sex   <- unique(pheno[, c("FID", "IID", "age", PCS), with = FALSE], by = c("FID", "IID"))
write_table(cov_with_sex, out("covars_with_sex.txt"))
write_table(cov_no_sex, out("covars_no_sex.txt"))

# Height-adjusted weight: sex, age, height, squared and cubed standardised height and log
# height, so that weight is adjusted for a flexible function of height, plus PC1-PC10 as for
# BMI (16 covariates). Height is standardised with the mean and SD over everyone with a height.
h_mu <- mean(pheno$height, na.rm = TRUE)
h_sd <- sd(pheno$height, na.rm = TRUE)
cov_weight <- unique(pheno[!is.na(height), .(FID, IID, sex, age, height,
                                             height_z2 = ((height - h_mu) / h_sd)^2,
                                             log_height = log(height),
                                             height_z3 = ((height - h_mu) / h_sd)^3)], by = c("FID", "IID"))
cov_weight <- merge(cov_weight, cov_with_sex[, c("FID", "IID", PCS), with = FALSE],
                    by = c("FID", "IID"), all.x = TRUE, sort = FALSE)
setcolorder(cov_weight, c("FID", "IID", "sex", "age", "height", "height_z2", "log_height", "height_z3", PCS))
write_table(cov_weight, out("covars_weight_height_model_pcs.txt"))

# Phenotype file for the GWAS of chromosomal sex (binary, male = 1).
write_table(cov_with_sex[, .(FID, IID, sex)], out("phenotypes_sex.txt"))

## 4. Keep files -------------------------------------------------------------------------
# A participant is kept when every covariate column is present and the phenotype(s) are
# present. Files list "FID FID", sorted as text. The analysed Stage 1 sample is this keep
# list further intersected with the REGENIE LOCO file and the BGEN sample file.
complete_ids <- function(cv) cv[complete.cases(cv), as.character(FID)]
pheno_ids <- function(cols) pheno[complete.cases(pheno[, ..cols]), as.character(FID)]
ids_with_sex <- complete_ids(cov_with_sex)
ids_no_sex   <- complete_ids(cov_no_sex)
keep_sets <- list(
  bmi     = intersect(ids_with_sex, pheno_ids("bmi")),
  log_bmi = intersect(ids_with_sex, pheno_ids("log_bmi")),
  weight  = intersect(ids_with_sex, pheno_ids("weight")),   # height terms are dropped as missing covariates later
  testosterone_series = intersect(ids_no_sex, pheno_ids(c("testosterone", paste0("testosterone_shrink_", c(0.2, 0.5, 0.7)))))
)
for (k in names(keep_sets)) {
  ids <- sort(unique(keep_sets[[k]]), method = "radix")
  write_keep(ids, ids, out("keep", paste0("keep_", k, ".txt")))
}
# Sex strata for the sex-stratified testosterone REGENIE (everyone with that sex code).
write_keep(cov_with_sex[sex == 1]$FID, cov_with_sex[sex == 1]$IID, out("keep", "keep_male.txt"))
write_keep(cov_with_sex[sex == 0]$FID, cov_with_sex[sex == 0]$IID, out("keep", "keep_female.txt"))

## 5. Antibody phenotypes (multiplex serology) -------------------------------------------
# Antigen levels are median fluorescence intensities (MFI) measured at up to two assessment
# instances. All antigens in the export are built the same way; the paper analyses JC virus
# VP1 (jcv_vp1) and Merkel cell polyomavirus VP1 (mcv_vp1). For each antigen:
#   MFI        = mean over the available instances
#   log_all    = log(MFI), all participants with MFI > 0           (Stage 1 phenotype)
#   seropos    = 1 if any instance reaches the assay cut-off, else 0 (cut-off 250 MFI for JCV/MCV VP1)
#   log_seropos = log(MFI) among seropositives only
# All antigens are kept in the build because REGENIE step 1 was fitted jointly to every
# all-sample log-MFI trait (02_regenie.sh), which fixes the sample used for the LOCO fit.
norm_name <- function(x) gsub("[[:space:]]+", " ", trimws(x))
slugify <- function(x) gsub("_+", "_", gsub("^_+|_+$", "", gsub("[^a-z0-9]+", "_", gsub("'", "", tolower(x)))))
# Antigens with an assay cut-off (antigen_id, column label in the export, cut-off in MFI).
antigens <- data.table(
  antigen_id = c("hsv1_mgg1", "hsv2_mgg2", "ebv_vca_p18", "ebv_ebna1", "ebv_zebra", "ebv_ead",
                 "cmv_pp150_nter", "cmv_pp52", "cmv_pp28", "hhv6_ie1a", "hhv6_ie1b", "hhv6_p101k",
                 "hhv7_u14", "vzv_ge_gi", "bkv_vp1", "jcv_vp1", "mcv_vp1",
                 "ctrach_momp_a", "ctrach_momp_d", "ctrach_tarpd_f1", "ctrach_tarpd_f2", "ctrach_porb", "ctrach_pgp3",
                 "hpylori_caga", "hpylori_vaca", "hpylori_omp", "hpylori_groel", "hpylori_catalase", "hpylori_urea",
                 "toxo_p22", "toxo_sag1"),
  pathogen_id = c("hsv1", "hsv2", rep("ebv", 4), rep("cmv", 3), rep("hhv6", 3), "hhv7", "vzv",
                  "bkv", "jcv", "mcv", rep("ctrach", 6), rep("hpylori", 6), rep("toxo", 2)),
  raw_mfi_base = c(
    "1gG antigen for Herpes Simplex virus-1", "2mgG unique antigen for Herpes Simplex virus-2",
    "VCA p18 antigen for Epstein-Barr Virus", "EBNA-1 antigen for Epstein-Barr Virus",
    "ZEBRA antigen for Epstein-Barr Virus", "EA-D antigen for Epstein-Barr Virus",
    "pp150 Nter antigen for Human Cytomegalovirus", "pp 52 antigen for Human Cytomegalovirus",
    "pp 28 antigen for Human Cytomegalovirus", "IE1A antigen for Human Herpesvirus-6",
    "IE1B antigen for Human Herpesvirus-6", "p101 k antigen for Human Herpesvirus-6",
    "U14 antigen for Human Herpesvirus-7", "gE / gI antigen for Varicella Zoster Virus",
    "BK VP1 antigen for Human Polyomavirus BKV", "JC VP1 antigen for Human Polyomavirus JCV",
    "MC VP1 antigen for Merkel Cell Polyomavirus", "momp A antigen for Chlamydia trachomatis",
    "momp D antigen for Chlamydia trachomatis", "tarp-D F1 antigen for Chlamydia trachomatis",
    "tarp-D F2 antigen for Chlamydia trachomatis", "PorB antigen for Chlamydia trachomatis",
    "pGP3 antigen for Chlamydia trachomatis", "CagA antigen for Helicobacter pylori",
    "VacA antigen for Helicobacter pylori", "OMP antigen for Helicobacter pylori",
    "GroEL antigen for Helicobacter pylori", "Catalase antigen for Helicobacter pylori",
    "UreA antigen for Helicobacter pylori", "p22 antigen for Toxoplasma gondii",
    "sag1 antigen for Toxoplasma gondii"),
  mfi_cutoff = c(150, 150, 250, 250, 100, 100, 100, 150, 200, 100, 100, 100, 100, 100,
                 250, 250, 250, 100, 100, 100, 100, 80, 200, 400, 100, 170, 80, 180, 130, 100, 160),
  has_cutoff = TRUE, order_group = 0L)
# Any other antigen column in the export is added as a log-MFI trait without a cut-off.
header <- names(fread(serology_file, nrows = 0))
inst <- data.table(col = grep(" \\| Instance [01]$", header, value = TRUE))
inst[, `:=`(base_norm = norm_name(sub(" \\| Instance [01]$", "", col)),
            instance = as.integer(sub("^.* \\| Instance ([01])$", "\\1", col)))]
other <- unique(inst[grepl("antigen", base_norm, ignore.case = TRUE) &
                     !grepl("assay QC|assay date", base_norm, ignore.case = TRUE) &
                     !base_norm %in% norm_name(antigens$raw_mfi_base), .(base_norm)])
if (nrow(other))
  antigens <- rbind(antigens, other[, .(antigen_id = paste0("raw_", slugify(base_norm)), pathogen_id = "raw_extra",
                                        raw_mfi_base = base_norm, mfi_cutoff = NA_real_, has_cutoff = FALSE,
                                        order_group = 1L)])
antigens[, raw_cols := lapply(raw_mfi_base, function(b) inst[base_norm == norm_name(b)][order(instance), col])]
setorder(antigens, order_group, pathogen_id, antigen_id)   # antigens with a cut-off first

ab <- fread(serology_file, select = unique(c("Participant ID", unlist(antigens$raw_cols))), na.strings = c("", "NA"))
setnames(ab, "Participant ID", "antibody_app_id")
for (col in setdiff(names(ab), "antibody_app_id")) set(ab, j = col, value = as.numeric(ab[[col]]))
ab[, antibody_app_id := as.integer(antibody_app_id)]
sl <- fread(serology_linker, select = c("ieu", "app"))
sl <- unique(sl[!is.na(ieu) & !is.na(app), .(FID = as.character(ieu), IID = as.character(ieu), antibody_app_id = as.integer(app))],
             by = "antibody_app_id")
# Covariates age, sex, PC1-PC10 (the run took them from an earlier export of the same fields).
ab_base <- unique(pheno[, c("FID", "IID", "age", "sex", PCS), with = FALSE][, `:=`(FID = as.character(FID), IID = as.character(IID))],
                  by = c("FID", "IID"))
ab <- merge(ab_base, merge(ab, sl, by = "antibody_app_id"), by = c("FID", "IID"))
mfi_cols <- intersect(setdiff(names(ab), c("FID", "IID", "antibody_app_id", "age", "sex", PCS)), unlist(antigens$raw_cols))
ab <- ab[rowSums(!is.na(ab[, ..mfi_cols])) > 0]          # anyone with at least one antigen measured

for (i in seq_len(nrow(antigens))) {
  id <- antigens$antigen_id[i]
  m <- as.matrix(ab[, intersect(antigens$raw_cols[[i]], names(ab)), with = FALSE])
  mfi <- rowMeans(m, na.rm = TRUE); mfi[is.nan(mfi)] <- NA_real_
  ab[, (paste0(id, "_mfi")) := mfi]
  ab[, (paste0(id, "_log_all")) := fifelse(!is.na(mfi) & mfi > 0, log(mfi), NA_real_)]
  if (antigens$has_cutoff[i]) {
    seen <- rowSums(!is.na(m)) > 0
    sero <- fifelse(seen, as.integer(rowSums(m >= antigens$mfi_cutoff[i], na.rm = TRUE) > 0), NA_integer_)
    ab[, (paste0(id, "_antigen_seropos")) := sero]
    ab[, (paste0(id, "_log_seropos")) := fifelse(sero == 1L & !is.na(mfi) & mfi > 0, log(mfi), NA_real_)]
  }
}
ab_covars <- unique(ab[, c("FID", "IID", "age", "sex", PCS), with = FALSE], by = c("FID", "IID"))
trait_cols <- grep("_(mfi|log_all|antigen_seropos|log_seropos)$", names(ab), value = TRUE)
write_table(ab[, c("FID", "IID", "antibody_app_id", trait_cols), with = FALSE], out("phenotypes_antibody_all.txt"))
write_table(ab_covars, out("covars_antibody_all.txt"))
complete_cov <- ab_covars[complete.cases(ab_covars[, c("age", "sex", PCS), with = FALSE]), .(FID, IID)]
write_keep(complete_cov$FID, complete_cov$IID, out("keep", "keep_antibody_all_complete.txt"))
# All-sample log-MFI traits with at least 100 values, in alphabetical order: this order is
# the phenotype order of the joint REGENIE step 1 and therefore the LOCO file index.
qt_all <- sort(paste0(antigens$antigen_id, "_log_all"), method = "radix")
qt_all <- qt_all[vapply(qt_all, function(t) sum(!is.na(ab[[t]])) >= 100, logical(1))]
writeLines(qt_all, out("antibody_qt_all_traits.txt"))
# Per-trait keep: complete trait and covariates (the Stage 1 keep for that trait).
for (t in c(qt_all, grep("_log_seropos$", names(ab), value = TRUE))) {
  m <- merge(ab[, c("FID", "IID", t), with = FALSE], ab_covars, by = c("FID", "IID"))
  k <- m[complete.cases(m[, c(t, "age", "sex", PCS), with = FALSE]), .(FID, IID)]
  write_keep(k$FID, k$IID, out("keep", paste0("keep_", t, ".txt")))
}

## 6. Trait-survey phenotypes ---------------------------------------------------------------
# Every survey trait is analysed in the frame of participants with complete age, sex and
# PC1-PC10, with the same covariates as BMI. A log column is written for positive traits;
# the scale analysed for each trait (raw or log) was chosen from the observed distribution
# by the Shapiro-Francia normality statistic W', computed here and saved as an aggregate table.
covariate_frame <- complete_ids(cov_with_sex)
shapiro_francia <- function(x) { n <- length(x); cor(sort(x), qnorm((seq_len(n) - 0.375) / (n + 0.25)))^2 }
moments <- list()

# 6a. Height, CRP and triglycerides from the core extract (positive values only).
d <- base[FID %in% covariate_frame]
survey_core <- d[, .(FID, IID)]
for (tr in c("height", "crp", "triglycerides")) {
  x <- as.numeric(d[[tr]]); ok <- is.finite(x) & x > 0
  v <- rep(NA_real_, nrow(d)); v[ok] <- x[ok]; survey_core[, (tr) := v]
  vl <- rep(NA_real_, nrow(d)); vl[ok] <- log(x[ok]); survey_core[, (paste0("log_", tr)) := vl]
  write_keep(d$FID[ok], d$IID[ok], out("keep", paste0("keep_", tr, ".txt")))
  moments[[tr]] <- data.table(trait = tr, n = sum(ok), W_raw = shapiro_francia(x[ok]), W_log = shapiro_francia(log(x[ok])))
}
write_table(survey_core, out("phenotypes_survey_core.txt"))

# 6b. Biochemistry, blood count, waist and spirometry panel. Field IDs (instance 0):
# 48 waist, 30760 HDL, 30780 LDL, 30870 triglycerides, 30790 lipoprotein(a), 30840 total
# bilirubin, 30730 GGT, 30620 ALT, 30750 HbA1c, 30880 urate, 30770 IGF-1, 30710 CRP,
# 30890 vitamin D, 30600 albumin, 30680 calcium, 30810 phosphate, 30860 total protein,
# 30150 eosinophil count, 30190 monocyte percentage, 30070 RDW, 30050 MCH, 30110 platelet
# distribution width, 20150 FEV1 (best), 20151 FVC (best); 6177/6153 medication (men/women).
# In the run, HbA1c (and the diabetes exclusion) came from this extract; the other panel traits'
# keep files were written by the same code from an export of the same fields from a second
# UK Biobank application (identical field list, linked to the same genetic identifiers).
panel_fields <- c(sex = 31, age = 21003, bmi = 21001, height = 50, weight = 21002, waist = 48, hdl = 30760,
                  ldl = 30780, triglycerides = 30870, lipoprotein_a = 30790, bilirubin_total = 30840, ggt = 30730,
                  alt = 30620, hba1c = 30750, urate = 30880, igf1 = 30770, crp = 30710, vitamin_d = 30890,
                  albumin = 30600, calcium = 30680, phosphate = 30810, total_protein = 30860,
                  eosinophil_count = 30150, monocyte_pct = 30190, rdw = 30070, mch = 30050,
                  platelet_dist_width = 30110, fev1_best = 20150, fvc_best = 20151)
hdr <- names(fread(ukb_extract, nrows = 0, sep = "\t"))
med_cols <- grep("^f\\.(6177|6153)\\.0\\.", hdr, value = TRUE)
diab_cols <- c("f.2443.0.0", "f.4041.0.0")
px <- fread(ukb_extract, select = c("f.eid", paste0("f.", panel_fields, ".0.0"), med_cols, diab_cols),
            sep = "\t", colClasses = "character")
setnames(px, c("f.eid", paste0("f.", panel_fields, ".0.0")), c("eid", names(panel_fields)))
# eid -> genetic ID, first occurrence of each eid; the linker CSV has quoted fields and CRLF endings.
lk <- fread(ukb_linker, colClasses = "character")
lk <- data.table(eid = lk[[4]], ieu = lk[[2]])[!duplicated(eid)]
px <- merge(px, lk, by = "eid", sort = FALSE)[!is.na(ieu) & ieu != ""]
# Self-reported cholesterol-lowering medication (code 1) in any instance-0 array of 6177/6153.
px[, statin_selfreport := as.integer(Reduce(`|`, lapply(med_cols, function(c) get(c) %in% "1")))]
panel <- px[, c("ieu", names(panel_fields), "statin_selfreport"), with = FALSE]
setnames(panel, "ieu", "FID")
for (c in names(panel_fields)) set(panel, j = c, value = suppressWarnings(as.numeric(fifelse(panel[[c]] == "", NA_character_, panel[[c]]))))
panel <- panel[FID %in% covariate_frame]
# FEV1/FVC is built as a ratio here so the construction is explicit.
panel[, fev1_fvc_ratio := fifelse(is.finite(fev1_best) & is.finite(fvc_best) & fvc_best > 0, fev1_best / fvc_best, NA_real_)]
survey_panel <- panel[, .(FID, IID = FID)]
for (tr in setdiff(names(panel), c("FID", "sex", "age", "statin_selfreport"))) {
  x <- panel[[tr]]; ok <- is.finite(x)
  if (sum(ok) < 5000) next
  v <- rep(NA_real_, nrow(panel)); v[ok] <- x[ok]; survey_panel[, (tr) := v]
  W_log <- NA_real_
  if (all(x[ok] > 0)) {
    vl <- rep(NA_real_, nrow(panel)); vl[ok] <- log(x[ok]); survey_panel[, (paste0("log_", tr)) := vl]
    W_log <- shapiro_francia(log(x[ok]))
  }
  write_keep(panel$FID[ok], panel$FID[ok], out("keep_panel", paste0("keep_", tr, ".txt")))
  moments[[tr]] <- data.table(trait = tr, n = sum(ok), W_raw = shapiro_francia(x[ok]), W_log = W_log)
}
write_table(survey_panel, out("phenotypes_survey_panel.txt"))
fwrite(rbindlist(moments), out("survey_scale_choice.csv"))   # aggregate only

# 6c. HbA1c without diabetes. Diabetes is defined from baseline self-report, never from HbA1c
# itself, because excluding on the outcome would truncate the distribution whose quantiles are
# estimated: field 2443 = 1 (doctor-diagnosed diabetes) unless field 4041 = 1 (gestational
# diabetes only), or insulin (code 3) in any instance-0 array of 6177/6153. "Do not know" and
# "prefer not to answer" are kept. The HbA1c column itself is the raw (unlogged) value.
px[, diabetes := (f.2443.0.0 %in% "1" & !(f.4041.0.0 %in% "1")) | Reduce(`|`, lapply(med_cols, function(c) get(c) %in% "3"))]
diab_ids <- px[diabetes == TRUE, ieu]
keep_hba1c <- fread(out("keep_panel", "keep_hba1c.txt"), header = FALSE, colClasses = "character")
keep_nodiab <- keep_hba1c[!V1 %in% diab_ids]
write_keep(keep_nodiab$V1, keep_nodiab$V2, out("keep_panel", "keep_hba1c_nodiab.txt"))

} else {

## 7. Samples defined after REGENIE step 1 ---------------------------------------------------
# The analysed Stage 1 sample is the keep list intersected with the participants in the
# REGENIE LOCO file (genotyped, passed step 1 QC) and in the imputed-data sample file.
loco_ids <- function(f) scan(text = readLines(f, n = 1L), what = character(), quiet = TRUE)[-1]   # header: FID_IID ...
in_loco <- function(fid, iid, ids) fid %in% ids | iid %in% ids | paste(fid, iid, sep = "_") %in% ids
bgen <- fread(bgen_sample, skip = 2L, header = FALSE, select = 1:2, colClasses = "character", col.names = c("FID", "IID"))
bgen_ids <- bgen$FID
rk <- function(f) fread(f, header = FALSE, colClasses = "character", col.names = c("FID", "IID"))
ph <- fread(out("phenotypes_paper.txt"), colClasses = list(character = c("FID", "IID")))

# 7a. BMI Stage 1 sample (leads-only Stage 1 runs in 03_stage1_quantile_gwas.sh).
k <- rk(out("keep", "keep_bmi.txt"))
k <- k[in_loco(FID, IID, loco_ids(loco_log_bmi)) & FID %in% bgen_ids]
write_keep(k$FID, k$IID, out("keep", "keep_bmi_stage1.txt"))

# 7b. Testosterone Stage 1 sample, for the GWAS of sex on exactly that sample.
cns <- fread(out("covars_no_sex.txt"), colClasses = list(character = c("FID", "IID")))
sx  <- fread(out("phenotypes_sex.txt"), colClasses = list(character = c("FID", "IID")))
k <- unique(rk(out("keep", "keep_testosterone_series.txt")))
k <- merge(k, ph[!is.na(testosterone), .(FID, IID)], by = c("FID", "IID"))
k <- merge(k, cns[complete.cases(cns), .(FID, IID)], by = c("FID", "IID"))
k <- k[paste(FID, IID, sep = "_") %in% loco_ids(loco_testosterone)]
k <- merge(k, bgen, by = c("FID", "IID"))
k <- merge(k, sx[!is.na(sex), .(FID, IID)], by = c("FID", "IID"))
write_keep(k$FID, k$IID, out("keep", "keep_testosterone_stage1.txt"))

# 7c. Antibody serostatus GWAS inputs on the Stage 1 antibody sample. Serostatus is the
# assay flag (any instance at or above 250 MFI); the seropositive-only trait is log MFI among
# seropositives. Values outside a trait's analysed sample are set to missing.
abph <- fread(out("phenotypes_antibody_all.txt"), colClasses = list(character = c("FID", "IID")))
abcv <- fread(out("covars_antibody_all.txt"), colClasses = list(character = c("FID", "IID")))
ka <- rk(out("keep", "keep_antibody_all_complete.txt"))
pred <- fread(antibody_pred_list, header = FALSE, col.names = c("trait", "loco"))
cols <- list(); cohorts <- list()
for (tr in c("jcv_vp1_log_all", "mcv_vp1_log_all")) {
  stem <- sub("_log_all$", "", tr)
  keep <- unique(rk(out("keep", paste0("keep_", tr, ".txt")))$FID)
  x <- abph[FID %in% keep & FID %in% abcv$FID]
  cvm <- abcv[match(x$FID, abcv$FID)]
  ok <- Reduce(`&`, lapply(setdiff(names(abcv), c("FID", "IID")), function(n) is.finite(as.numeric(cvm[[n]]))))
  x <- x[ok & is.finite(as.numeric(get(tr)))]
  x <- x[in_loco(FID, IID, loco_ids(pred[trait == tr]$loco)) & FID %in% bgen_ids & FID %in% ka$FID]
  s <- as.integer(x[[paste0(stem, "_antigen_seropos")]])
  d <- data.table(FID = x$FID, IID = x$IID, s = s, lp = ifelse(s == 1, as.numeric(x[[tr]]), NA_real_))
  setnames(d, c("s", "lp"), c(paste0(stem, "_seropos"), paste0(stem, "_log_seropos")))
  cols[[tr]] <- d; cohorts[[tr]] <- x[, .(FID, IID)]
  write_keep(x$FID[s == 1], x$IID[s == 1], out("keep", paste0("keep_", tr, "_seropos.txt")))
}
write_table(Reduce(function(a, b) merge(a, b, by = c("FID", "IID"), all = TRUE), cols), out("phenotypes_serostatus.txt"))
ka2 <- unique(rbindlist(cohorts))
write_keep(ka2$FID, ka2$IID, out("keep", "keep_antibody_analysed.txt"))
}
