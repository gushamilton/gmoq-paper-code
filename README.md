<p align="center"><img src="docs/gmoq_logo.png" alt="G-MoQ" width="380"></p>

# Code for "Dissecting complex genetic effects and phenotype scale with the Genomic Method of Quantiles"

Hamilton, Butler-Laporte, Pires Hartwig, Hemani, Tilling, Davey Smith, Nivard.

This repository contains the analysis scripts, aggregate results, figures and supplementary tables for the paper. The method itself—the Genomic Method of Quantiles (G-MoQ)—is a separate software package: <https://github.com/gushamilton/GMoQ>. It has two stages:

- **Stage 1** (`gmoq-stage1`, Python/C++) runs an unconditional quantile GWAS. Each variant's effect is estimated at 45 quantiles of the phenotype ($\tau = 0.02, \dots, 0.98$) by recentred influence function (RIF) regression, with a pooled covariance of the 45 effects.
- **Stage 2** (`gmoq`, R) converts each variant's 45 quantile effects into effects on the parameters of a user-specified model (such as a two-component mixture or a left-censored Tobit model), or summarises them without a model via moments and scale-invariant whitened principal components.

The analyses use UK Biobank data (application 81499), which are available to approved researchers from UK Biobank. This repository contains only the code and aggregate summary statistics needed to reproduce every figure and table (`results/`). No individual-level or participant-identifiable data are included.

## Rebuilding the figures and tables

Every main figure, Extended Data figure, supplementary figure and supplementary table is rebuilt directly from `results/` in under a minute:

```bash
make all
```

This requires R 4.5 with `data.table`, `ggplot2`, `patchwork`, `ragg` and `yaml`. Outputs are written to `paper_assets/`: PNG and PDF artwork, source data CSVs and legends.

## Where is the code for each method?

| If you want to… | Look here |
|---|---|
| Run Stage 1 quantile GWAS or Stage 2 parametric models | <https://github.com/gushamilton/GMoQ> (examples: `analysis/02_mixtures/01_fit_mixtures.R`, `02_genomewide_stage2.R`) |
| Run the **scale-invariant whitened PCA**, parallel analysis and per-variant **PC2 scores** | `R/pca.R` and `analysis/06_traits/01_moments_and_principal_components.R` |
| Convert Stage 1 profiles across phenotype transformations (raw, log, Box–Cox) | `R/scale.R` and `analysis/06_traits/04_change_of_scale.R` |
| Compute L-moments and between-lead heterogeneity | `R/moment_decomposition.R` and `R/heterogeneity.R` |
| Reproduce every figure (Figs. 1–7, Extended Data Fig. 1, S1–S15) and table (S1–S8) | `make all` (`figures/` and `tables/`) |

## How the repository follows the paper

| Paper section | What is done | Analysis scripts | Figure / table scripts |
|---|---|---|---|
| Methods: data and GWAS | Phenotypes, REGENIE GWAS, Stage 1 quantile GWAS, LD score regression, lead variants, BMI polygenic score | `analysis/00_gwas/` | |
| Methods: shared inputs | Baseline quantiles, lead-variant sets, lead quantile profiles, Stage 1 mean against REGENIE | `analysis/01_data/` | Supplementary Fig. S12, Table S7 |
| Simulations | Simulation studies and speed benchmarks | `analysis/08_simulations/` | Supplementary Figs. S1–S6, Tables S1–S2 |
| Testosterone: sex-blind mixture & censoring | Conceptual workflow, mixture fit and Stage 2 weights, genome-wide Stage 2, sex-stratified concordance, speed benchmarks, limit-of-detection Tobit censoring | `analysis/02_mixtures/`, `analysis/03_testosterone/` | Figs. 1, 2; Extended Data Fig. 1; Supplementary Figs. S5, S7, S8; Table S3 |
| Antibodies: membership versus level | Mixture decomposition of JC and Merkel cell polyomavirus IgG associations; thresholded serostatus GWAS | `analysis/04_antibodies/` | Fig. 3; Supplementary Fig. S9; Table S4 |
| BMI: competing models and polygenic background | Model profiles and angles; per-lead model scores; leave-one-chromosome-out polygenic-score interactions | `analysis/05_bmi/` | Fig. 4; Supplementary Figs. S13–S15; Table S6 |
| Scale-invariant dimensionality across 21 traits | Whitened PCA, parallel analysis, L-moments, between-lead heterogeneity, change of scale, SIQReg comparison, and biological dissection of 2D traits via scale-invariant PC2 scores | `analysis/06_traits/` | Figs. 5, 6, 7; Supplementary Fig. S10; Table S5 |
| HbA1c and diabetes | Third scale-invariant dimension, skew beyond shift and spread, and blind diabetes mixture | `analysis/07_hba1c/` | Figs. 5d, 6c; Supplementary Fig. S11 |
| Software | Versions used | `analysis/09_software/` | Supplementary Table S8 |

## Layout

| Path | Contents |
|---|---|
| `config.yaml` | Analysis settings: quantile grid, estimator choices, trait list, sample sizes, lead-variant rules |
| `R/` | Shared helpers (`source("R/load.R")` loads them all) |
| `analysis/` | Cluster analysis scripts, one folder per section of the paper |
| `results/<section>/<analysis>/` | Aggregate summary outputs read by the figure and table scripts |
| `figures/`, `tables/` | One self-contained script per figure or supplementary table |
| `paper_assets/` | Rendered figures and tables, with source data and legends |

## Running the analyses on your own data

1. Install G-MoQ from <https://github.com/gushamilton/GMoQ> (`pip install ./stage1` and `R CMD INSTALL stage2`).
2. Prepare phenotypes and covariates (`analysis/00_gwas/`), update the paths in `config.yaml`, and run the section folders in `analysis/` in order.
3. Run `make all` to regenerate all figures and tables in `paper_assets/`.

## Licence

MIT for the code. Aggregate results derive from UK Biobank (application 81499; the trait survey also used application 43612) and contain no individual-level data.
