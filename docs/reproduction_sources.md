# Reproduction sources for BMI model comparison and biological correlates

| Display or result | Figure/table producer | Upstream producer | Released aggregate inputs |
|---|---|---|---|
| Figure 4 | `figures/fig4_bmi_models.R` | `analysis/05_bmi/01_model_profiles.R`; lead-profile extraction | `results/05_bmi/model_profiles/`; `results/01_data/lead_profiles/` |
| Figure 4 model-separation claim | Figure 4 legend | `analysis/05_bmi/03_model_separation.R` | `results/05_bmi/model_separation/`; raw/log lead profiles, pooled covariances and model weights |
| Figure 7 ALT and LDL | `figures/fig7_hidden_biology.R` | `analysis/06_traits/biological_correlates/run_hidden_biology.py` and `.sh` | `results/06_traits/hidden_biology/alt_*`, `ldl_*`; primary lead moments |
| Figure 7 urate and lung function; CRP and bilirubin in Table 1 | Figure 7 and Table 1 producers | `analysis/06_traits/biological_correlates/run_more_traits.py` and `.sh` | `results/06_traits/hidden_biology/urate_*`, `fev_*`, `crp_*`, `bilirubin_*`; primary lead moments |
| Main Table 1 | `tables/table1_biological_correlates.R` | The biological-correlate producers above; primary moment/PCA analysis; HbA1c diabetes analysis | `results/06_traits/moments/`; `results/06_traits/hidden_biology/`; `results/07_hba1c/skew_diabetes/` |

The Figure 4 producer uses the manuscript's original analytic construction. The
angle matrix includes the infinitesimal mean-matched log-variance profile, with
the chromosome-1 raw-BMI pooled covariance defining its whitened metric. The
source CSV exports the six calculated pairwise angles in both metrics and both
reference curves actually plotted in panel c. The legacy pairwise-angle CSV is
not substituted for these calculations.

The BMI comparison was recovered from the original aggregate-only producer,
retaining the three equal-rank model pairs, covariance, mean matching, simulation
seed and draw order. Its saved full simulation used 20,000 draws. The port was
validated with a 20-draw smoke run to a temporary output directory; residual
scores matched the original model scores to less than 1e-12 relative error.
The saved full-calibration per-lead results were reused. No participant analysis
or full simulation rerun was performed to prepare these sources.

The original biological-correlate Python source files were recovered from their
cluster input directories. Their numerical function syntax trees are unchanged
(11 functions in the ALT/LDL producer; 12 in the other-trait producer). Changes
supply paths from environment variables, add reproduction instructions, and
correct descriptions of the historical subgroup procedure. The original source
SHA-256 values are recorded below; no restricted data, identifiers or logs were
copied. The original scripts did not freeze a Python dependency environment.
NumPy, pandas and SciPy are required; exact historical package versions were not
recorded in those source files.

The primary Table 1 ratios are taken directly from the Figure 6b
`pca_parallel.csv` rows: observed PC2 share divided by the one-shape-null median.
The subgroup nulls in the original auxiliary producers differ from that primary
procedure, as documented in
`analysis/06_traits/biological_correlates/README.md`. Table 1 retains the
interaction and clinical correlations and the LDL profile/asymmetry comparison,
and makes no claim about a change in primary dimension after restriction.

Table 1 writes a display CSV, caption and full-precision numeric source CSV.
It keeps the score definition, lead denominator and nominal two-sided P value
for each correlation. The whole-cohort and no-diabetes HbA1c asymmetry counts
use log-scale lead panels; the no-diabetes component ratio and its PC2
correlation use the raw-scale panel.

| Original source | SHA-256 |
|---|---|
| `run_hidden_biology.py` | `d0cba2545872ee7ebdc74d84c773d899816d7d76f14320d73e30062214ee4923` |
| `run_more_traits.py` | `9ebd6054249a2c704a8a8be395fd6242f43af85771cb18a6d022be738cf1bb17` |
