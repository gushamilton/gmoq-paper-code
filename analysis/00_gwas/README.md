# 00_gwas: from UK Biobank to the GWAS inputs

These scripts show how the genome-wide inputs of the paper were made from UK Biobank. They need
individual-level UK Biobank data and a compute cluster, so they are not run from this repository;
every path is a placeholder (`/path/to/...`) to replace with your own. Their outputs are
participant-level or genome-wide files that are not distributed. The analyses in `analysis/01_data`
onwards start from aggregate summaries of these outputs in `results/`.

Software: REGENIE v4.1.2, PLINK 2.00a6.8LM, gmoq-stage1 (G-MoQ Stage 1, commit in
`config.yaml`), LDSC 2.0.0, bgenix 1.1.7, R 4.5.1.
Genotypes: UK Biobank directly genotyped array (REGENIE step 1) and imputed BGEN (everything else);
HapMap3 variants (1,200,240 after intersection). Quantile levels: 45, `linspace(0.02, 0.98, 45)`.

## Steps, in order

| script | does | produces |
|---|---|---|
| `01_build_phenotypes.R` | links UK Biobank fields to the genetic IDs; builds BMI, log BMI, height-adjusted weight (16 covariates), log testosterone and its shrunk arms (F = 0.2, 0.5, 0.7), JCV/MCV VP1 antibody levels and serostatus, survey traits and HbA1c without diabetes | phenotype, covariate and keep files |
| `02_regenie.sh` | REGENIE step 1 (LOCO predictions) and step 2 (mean GWAS): BMI, log BMI, weight, sex-stratified testosterone, antibodies (all-sample, serostatus, seropositive-only), GWAS of sex on the testosterone Stage 1 sample, survey step 1 | LOCO prediction lists; `.regenie` files |
| `01_build_phenotypes.R --after-step1` | samples that depend on the LOCO files: BMI and testosterone Stage 1 samples, antibody serostatus inputs | keep files, serostatus phenotype file |
| `03_stage1_quantile_gwas.sh` | Stage 1 quantile GWAS per trait and chromosome, with the flags of each trait group; leads-only runs for the polygenic-score analyses | `scan.stage1.tsv.gz`, `scan.sigma0.tsv.gz`, `rif_nuisance.json` per chromosome |
| `04_ldsc.sh` | LDSC munging, SNP heritability and genetic correlation | LDSC logs |
| `05_lead_variants.sh` | plink2 clumping (P < 5e-8, secondary P < 1e-3, r2 < 0.01, 10 Mb, 1000 Genomes EUR) | lead variant lists |
| `06_bmi_prs_clump_threshold.sh` | BMI polygenic score by PLINK clumping and thresholding (REGENIE BMI GWAS, r2 < 0.1 within 250 kb, 1000 Genomes EUR, index P < 1e-5), one partial score per chromosome | partial scores per chromosome |

Stage 1 settings by trait group (all: 1000 jackknife blocks, batch size 100, seed 42):

| traits | covariates | LOCO | covariance output | variants |
|---|---|---|---|---|
| BMI, log BMI | age, sex, PC1-10 | covariate | full | HapMap3 |
| testosterone (4 arms) | age, PC1-10 | covariate | full | HapMap3 |
| height-adjusted weight | sex, age, height, height_z², log height, height_z³, PC1-10 | covariate | pooled | HapMap3 |
| JCV/MCV VP1 log MFI | age, sex, PC1-10 | covariate (regression scale) | pooled-only | HapMap3 with in-sample MAF > 0.05 |
| survey traits | age, sex, PC1-10 | covariate for height, log CRP, log triglycerides; none otherwise | pooled | HapMap3 |

## UK Biobank fields used

Instance 0 (baseline) throughout.

- Core: 31 sex, 21003 age, 50 standing height, 21002 weight, 21001 BMI, 30850 testosterone,
  22009 genetic principal components 1-10.
- Survey: 48 waist, 30760 HDL cholesterol, 30780 LDL cholesterol, 30870 triglycerides,
  30790 lipoprotein(a), 30840 total bilirubin, 30730 GGT, 30620 ALT, 30750 HbA1c, 30880 urate,
  30770 IGF-1, 30710 CRP, 30890 vitamin D, 30600 albumin, 30680 calcium, 30810 phosphate,
  30860 total protein, 30150 eosinophil count, 30190 monocyte percentage, 30070 red cell
  distribution width, 30050 mean corpuscular haemoglobin, 30110 platelet distribution width,
  20150 FEV1 (best measure), 20151 FVC (best measure).
- Medication and diabetes: 6177 and 6153 (medication for cholesterol, blood pressure, diabetes or
  hormones; code 1 cholesterol lowering, code 3 insulin), 2443 doctor-diagnosed diabetes,
  4041 gestational diabetes only.
- Serology: the multiplex serology antigen MFI columns of the UK Biobank infectious-disease
  serology export (instances 0 and 1), read by antigen name; JC virus VP1 and Merkel cell
  polyomavirus VP1 are analysed.
