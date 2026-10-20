<p align="center"><img src="docs/gmoq_logo.png" alt="G-MoQ" width="380"></p>

# Code for "Flexible, functional GWAS of complex effects using fast unconditional quantile regression"

Hamilton, Pires Hartwig, Hemani, Tilling, Davey Smith, Nivard.

This repository contains the analyses, figures and tables of the paper. The method itself,
the Genomic Method of Quantiles (G-MoQ), is a separate software package:
<https://github.com/gushamilton/GMoQ>. It has two stages:

- **Stage 1** (`gmoq-stage1`, Python/C++) runs an unconditional quantile GWAS. Each variant's
  effect is estimated at 45 quantiles of the phenotype (τ = 0.02, …, 0.98) by recentred
  influence function (RIF) regression, with a pooled covariance of the 45 effects.
- **Stage 2** (`gmoq`, R) converts each variant's 45 quantile effects into effects on
  the parameters of a model, such as the two means and the membership of a mixture. It can
  also summarise them without a model, by moments or principal components.

The analyses use UK Biobank data (application 81499), which are available to approved
researchers from UK Biobank. This repository contains the code and the summary results
needed to reproduce the figures and tables (`results/`: summary statistics, fitted model
parameters and aggregate counts).

## Rebuilding the figures and tables

Every figure and supplementary table is rebuilt from `results/` in under a minute:

```bash
make all
```

This needs R 4.5 with data.table, ggplot2, patchwork, ragg and yaml. The outputs go to
`paper_assets/`: a PNG and PDF, the source data and the legend for each item.

## How the repository follows the paper

| Paper | What is done | Analysis scripts | Figure / table scripts |
|---|---|---|---|
| Methods: data and GWAS | Phenotypes, REGENIE GWAS, Stage 1 quantile GWAS, LD score regression, lead variants, BMI polygenic score | `analysis/00_gwas/` | |
| Methods: shared inputs | Baseline quantiles, lead-variant sets, lead quantile profiles, Stage 1 mean against REGENIE | `analysis/01_data/` | Supplementary Fig. S12, Table S7 |
| Simulations | The simulation studies (drivers for studies 1–7) and the import of their summaries | `analysis/08_simulations/` | Supplementary Figs. S1–S6, Tables S1–S2 |
| Testosterone: a sex-blind GWAS read through a two-component mixture | Mixture fit and Stage 2 weights (model and hybrid density); genome-wide Stage 2; comparison with sex-stratified GWAS and a GWAS of sex; exemplar loci; speed | `analysis/02_mixtures/`, `analysis/03_testosterone/` | Figs. 1, 2; Supplementary Figs. S5, S7, S8; Table S3; Supplementary Note 2 |
| Antibodies: membership versus level | Mixture decomposition of JC and Merkel cell polyomavirus antibody associations; thresholded serostatus GWAS; FUT2 threshold analysis | `analysis/04_antibodies/` | Fig. 3; Supplementary Fig. S9; Table S4; Supplementary Note 3 |
| BMI: competing models, polygenic background | Model profiles and angles; per-lead model scores; leave-one-chromosome-out polygenic-score interactions | `analysis/05_bmi/` | Figs. 4, 6; Supplementary Figs. S13, S14; Table S6 |
| Patterns across 21 trait analyses | L-moments, lack of fit, principal components, between-lead heterogeneity, change of scale, comparison with SIQReg | `analysis/06_traits/` | Fig. 5a,b; Supplementary Fig. S10; Table S5 |
| HbA1c and diabetes | Skew beyond shift and spread, and a mixture fitted without knowledge of diabetes | `analysis/07_hba1c/` | Fig. 5c; Supplementary Fig. S11 |
| Software | Versions used | `analysis/09_software/` | Supplementary Table S8 |

Within each folder the scripts are numbered in the order they run. Each one starts with a
short header that says what it does, which figure or table it supports, what it reads and
writes, and where it runs. Scripts marked "locally from results/" run on a laptop. Scripts
marked "on a cluster with UK Biobank data" need your own data.

## Layout

| Path | Contents |
|---|---|
| `config.yaml` | Every setting: paths to the UK Biobank-derived inputs, the quantile grid, estimator choices, trait list and sample sizes, lead-variant rules |
| `R/` | Shared functions: L-moment functionals and lack of fit, whitening and principal components, between-lead heterogeneity, change of scale, mixtures, reading GWAS files, reading and writing results. `source("R/load.R")` loads them all |
| `analysis/` | The analyses, one folder per part of the paper (table above) |
| `results/<section>/<analysis>/` | Aggregate outputs of the analyses, read by the figures and tables |
| `figures/`, `tables/` | One script per figure or supplementary table, named after it |
| `paper_assets/` | Rendered figures and tables, with source data and legends |

## Running the analyses with your own data

1. Install G-MoQ, `gmoq` 0.7.0 or later (Python package and command `gmoq-stage1`; R package
   `gmoq`). The paper's analyses ran on development versions of the same
   code, before the packages took their current name; `gmoq` 0.7.0 has the same estimators.
   Some defaults have changed since (for example, mixture membership weights now default to
   the hybrid density), so set the options that `config.yaml` records.
2. Build the GWAS inputs with the scripts in `analysis/00_gwas/`; its README lists the steps
   and the UK Biobank fields.
3. Edit `config.yaml` so that its paths point to your files. The paths in it are those of
   the authors' cluster.
4. Run `analysis/01_data/` and then the section folders in order. Outputs are written to
   `results/`, or to `$GMOQ_OUT` if it is set, which is convenient for comparing a rerun with
   the published results.
5. Run `make all`.

## Licence

MIT for the code. The results derive from UK Biobank (application 81499; the trait survey
also used application 43612) and are shared as aggregate statistics only.
