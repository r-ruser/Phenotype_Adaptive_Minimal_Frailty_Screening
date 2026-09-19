# Phenotype-Adaptive Minimal Frailty Screening

Reproducible R workflow for recovering harmonized same-wave data, discovering clinical-social phenotypes, deriving phenotype-adaptive 4/5-item frailty screeners, performing leave-one-cohort-out validation, and exporting Nature-style figures.

## Analysis sequence

1. `scripts/01_recover_same_wave_data.R`
2. `scripts/02_lca_transportability.R`
3. `scripts/03_irt_dif.R`
4. `scripts/04_screening_validation.R`
5. `scripts/05_nature_figures_and_report.R`

Run the complete workflow with:

```powershell
Rscript scripts/run_all.R
```

Required R packages: `haven`, `dplyr`, `tidyr`, `poLCA`, `ltm`, `pROC`,
`ggplot2`, `patchwork`, `scales`, `svglite`, and `ragg`.

Raw harmonized ageing-study files remain outside this repository. Scripts resolve them from the configured local data root. Generated participant-level data are not committed.

## Reproducibility contract

- Baseline and follow-up variables are extracted from explicitly paired waves within each cohort.
- Survey negative/special codes are converted to missing values before derivation.
- Phenotype discovery uses baseline context only.
- Screening outcomes and Paper 1 treatment-benefit estimates do not enter phenotype discovery or item selection.
- Leave-one-cohort-out validation locks item selection, model fitting, and screening thresholds in training cohorts.
- Figure exports include SVG, PDF, TIFF, PNG, and panel-level source data.

## Generated deliverables

The workflow overwrites the analysis-ready RDS files and tabular results in the
parent analysis directory. Manuscript-facing outputs include
`Paper2_Analysis_Summary.md`, protocol-aligned `output/TABLE1_*.csv` through
`output/TABLE6*.csv`, and the complete `output/nature_figures/` bundle. Large
data and result files are intentionally excluded from Git.
