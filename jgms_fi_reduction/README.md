# Five-cohort FI item-reduction methods package

This package supplies the actual coding tables, aggregate results, and configurable R analysis code used for the five-cohort study. The source repository is [Phenotype_Adaptive_Minimal_Frailty_Screening](https://github.com/r-ruser/Phenotype_Adaptive_Minimal_Frailty_Screening).

The primary score uses a universal training-only ranking. The six-item form is an exploratory, practical candidate length identified by describing the held-out length curves. The package records the distinction between training-only item ranking and descriptive length selection.

## Package contents

- `mappings/complete_FI_encoding.csv`: all 222 included cohort-item mappings, covering 46 HRS, 44 ELSA, 50 SHARE, 39 CHARLS, and 43 MHAS deficits; source field names, labels, response values, and exact 0–1 coding rules.
- `mappings/reference_item_mapping_all_audited.csv`: both retained and screened-out reference items, with availability and inclusion flags.
- `mappings/candidate16_encoding.csv`: the full 16-item candidate bank in all five cohorts.
- `mappings/FI6_reader_scoring_table.csv` and `code/score_FI6.R`: a reader-facing six-item scoring specification and usable scoring functions.
- `mappings/incident_ADL_definition.csv`: the exact baseline/follow-up aggregate variables and task composition used in the primary incident ADL analysis.
- `mappings/ADL_codebook_evidence.csv`: brief paraphrases with local official codebook names and PDF page numbers.
- `results/LCA_*`: every candidate class-number fit, selected classes, complete-context fitting denominators, class sizes, and conditional profiles.
- `results/actual_wave_years.csv`: observed baseline/follow-up years and interview intervals.
- `results/`: aggregate rankings, thresholds, length curves, paired comparisons, overlap sensitivity, and verification records.
- `documentation/supplementary_methods.md`: detailed methods based on inspected code and actual data.
- `documentation/input_specification.md`: source data, intermediate object schemas, data access, and processing order.
- `documentation/JGMS_position_and_historical_data.md`: article positioning and a cover-letter paragraph explaining historical waves.

Participant records, raw survey files, individual-level scores, and participant identifiers remain in the user's private data environment. Each licensed researcher obtains the source files through the cohort's access process.

## Configure and run

The source data root contains the `HRS_USA`, `ELSA_UK`, `SHARE_Europe`, `CHARLS_China`, and `MHAS_Mexico` folders specified in the input document. `MS5_RUNTIME_DIR` is a private output folder outside the public repository. Set the package path and the two locations in PowerShell:

```powershell
$env:MS5_PACKAGE_DIR = 'C:/research/methods_package'
$env:MS5_DATA_ROOT = 'E:/licensed_cohort_data'
$env:MS5_RUNTIME_DIR = 'E:/private_analysis_runtime'
& 'D:/R-4.4.3/bin/x64/Rscript.exe' "$env:MS5_PACKAGE_DIR/code/01_recover_baseline.R"
& 'D:/R-4.4.3/bin/x64/Rscript.exe' "$env:MS5_PACKAGE_DIR/code/02_prepare_complete_FI.R"
& 'D:/R-4.4.3/bin/x64/Rscript.exe' "$env:MS5_PACKAGE_DIR/code/03_candidate_bank.R"
& 'D:/R-4.4.3/bin/x64/Rscript.exe' "$env:MS5_PACKAGE_DIR/code/04_training_phenotypes_SHAP.R"
& 'D:/R-4.4.3/bin/x64/Rscript.exe' "$env:MS5_PACKAGE_DIR/code/05_universal_SHAP.R"
& 'D:/R-4.4.3/bin/x64/Rscript.exe' "$env:MS5_PACKAGE_DIR/code/06_locked_scores.R"
& 'D:/R-4.4.3/bin/x64/Rscript.exe' "$env:MS5_PACKAGE_DIR/code/07_performance.R"
& 'D:/R-4.4.3/bin/x64/Rscript.exe' "$env:MS5_PACKAGE_DIR/code/08_retention.R"
& 'D:/R-4.4.3/bin/x64/Rscript.exe' "$env:MS5_PACKAGE_DIR/code/09_paired_comparisons.R"
& 'D:/R-4.4.3/bin/x64/Rscript.exe' "$env:MS5_PACKAGE_DIR/code/10_overlap_sensitivity.R"
```

Required R packages are `haven`, `dplyr`, `pROC`, `xgboost`, and `glmnet`; figure extensions additionally use `ggplot2`, `patchwork`, `svglite`, `ragg`, `xml2`, and `pdftools`. R 4.4.3 was used for the verification in this package. The fixed model settings and random seeds are in the methods and code; the attached session record specifies the available runtime. The two preparation stages were checked against the frozen eligible cohort, continuous FI, and primary ADL inputs using an in-memory audit. The training/model chain is supplied as executable code with its original aggregate outputs; a complete model refit belongs to the licensed source-data environment.


## Revision extensions

After stages 01–10, run `11_continuous_agreement.R`, `12_six_item_benchmarks.R`, and `13_draw_benchmark_figure.R` for continuous agreement, paired retention, and selection benchmarks. Run `14_SHARE_country_source.R` before `15_attrition_IPW.R`; then run `16_shared_five_ADL_sensitivity.R`. They read private intermediate objects and create outputs under `MS5_RUNTIME_DIR/extensions`. Aggregate results from the executed original analyses are included in `results/approximation`, `results/benchmarks`, and `results/attrition`. The path-configurable copies were syntax checked; a complete portable model-chain refit was not repeated for this revision.

The first two preparation stages were numerically verified against frozen analysis inputs. Supplementary intervals condition on the locked selected items; the six-item length remains an exploratory candidate. Training-based linear mapping is supplementary scale calibration. Random-set distributions describe variation across combinations. Outcome-availability weighting targets explicitly alive participants in countries with observed follow-up FI; recorded deaths and unresolved vital statuses are separate populations.

No participant data are included in this public package. Intermediate RDS files and the participant risk-set CSV created during local analysis belong to the private runtime directory.

## Additional result documentation

Figures S1 and 1 are reproduced with `17_draw_confusion_matrix.R` and `18_draw_framework.R`, using the supplied aggregate classification table and explicit node/edge specification. The framework script also uses `systemfonts`. These drawing scripts and the benchmark drawing script were executed in a separate runtime during revision packaging.
