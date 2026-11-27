# Biological correlates of component scores

`run_hidden_biology.py` and `run_more_traits.py` are the original numerical
producers of the ALT, LDL, urate, lung-function, CRP and bilirubin summaries.
Their restricted input paths are supplied through environment variables; their
numerical functions have been preserved. The accompanying shell scripts preserve
the indexed BGEN extraction, allele-aligned dosage export and binary-trait GLMs.
Only source code and aggregate variant or trait summaries are included here.

`tables/table1_biological_correlates.R` reproduces Main Table 1 from these
summaries and the released primary moment/parallel-analysis and HbA1c results.
It also writes full-precision numeric source values. The scores are different:
HbA1c uses unit-length whitened PC2 directions; Figure 7a-c, CRP and bilirubin use
unnormalised whitened PC2 projections (`whitened_norm * pc2_score`). The
correlation summaries originally emitted by these Python programs use their input
`pc2_score`, which is unit-normalised. Use the Table 1 producer for the
unnormalised correlations reported in the manuscript.

## Historical subgroup procedure

The subgroup calculations in these original producers reduce the rank of the
unit-normalised whitened profiles and then restore row lengths when fitting the
parallel-analysis null. The primary procedure in `R/pca.R` fits that null to
unnormalised whitened profiles before unit normalisation. These operations are
not equivalent. The original producers also use NumPy's default quantile rule
(type 7), whereas the primary Stage 1 analysis used type 8. The saved subgroup
component ratios and counts are retained as provenance, and are not used to
claim that primary dimensionality changed after a restriction. Table 1 takes all
PC2 ratios from the primary Figure 6b source. The original LDL profile comparison
and the six-to-zero asymmetry count remain descriptive results across the same
230 leads.

`masld_any` is a historical column name for a composite recorded-liver-disease
outcome, including fatty liver and fibrosis/cirrhosis codes. It does not diagnose
occult MASLD. `statin` is a historical column name for reported cholesterol-lowering
medication (fields 6177/6153, code 1); the restricted sample has no reported cholesterol-lowering medication and
includes participants with missing medication status. Bilirubin `dom_dev_z` tests the quadratic dosage term `(g - 1)^2`
adjusted for additive dosage, using imputed dosages.

## Restricted-data reproduction

These programs require an authorised UK Biobank environment and a BP compute
node. They do not run when rebuilding figures/tables. Python requires NumPy,
pandas and SciPy; shell orchestration uses bgenix and PLINK 2. The original
PLINK module was `apps/plink2/2.00a68LM`. Supply the original application-81499
extract files `data.51913.csv`, `data.669894.csv` (the latter only for the ALT/LDL
producer), and `linker.81499.csv`, together with the matching covariate table,
cohort keep lists and BGEN sample. The original extract field definitions and
inclusion rules are visible in `run_prep()` in each producer. Identically named
extract files from a different release need not have identical contents.

Export these environment variables in the allocated compute shell:

- `GMOQ_UKB_DIR`: directory holding the named restricted application extracts.
- `GMOQ_COVARIATES`: matching age, sex and ten-PC covariate table.
- `GMOQ_BGEN_SAMPLE`: imputed-genotype sample file.
- `GMOQ_BGEN_DIR`: chromosome BGEN directory, with `data.chrNN.bgen` and indexes.
- `GMOQ_PRIVATE_ROOT`: private output directory outside the public repository.
- `GMOQ_KEEP_ALT` and `GMOQ_KEEP_LDL`: matching ALT/LDL keep files for the first producer.
- `GMOQ_KEEP_DIR`: directory with `keep_fev1_fvc_ratio.txt`, `keep_urate.txt`,
  `keep_bilirubin_total.txt`, and `keep_crp.txt` for the second producer.
- `BGENIX`: bgenix executable; `PYTHON`: optional Python executable.
- `PLINK_MEMORY_MB`: memory available to PLINK within the allocation.

Select a separate private output root for each producer. Prepare its aggregate
lead inputs, then run its wrapper from the repository root:

```bash
python3 analysis/06_traits/biological_correlates/prepare_lead_inputs.py \
  --analysis hidden --root "$GMOQ_PRIVATE_ROOT"
bash analysis/06_traits/biological_correlates/run_hidden_biology.sh
```

For the other producer, use `--analysis more` and `run_more_traits.sh`. Allocate
resources after a representative smoke test on the intended cohort and hardware;
the wrappers request no allocation themselves. The genotype matrices, phenotype
rows and 45-column RIF matrix scale with the cohort size and lead count. All
participant-level files stay under `GMOQ_PRIVATE_ROOT`. Only reviewed files in
its `aggregate/` directory may be exported. The released copies of those
aggregates are in `results/06_traits/hidden_biology/`.
