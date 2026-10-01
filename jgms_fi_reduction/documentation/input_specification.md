# Licensed inputs and processing specification

## Source files

Set `MS5_DATA_ROOT` to the parent of the five cohort folders. The expected relative paths are:

| Cohort | Expected file | Release evidence | Join key |
|---|---|---|---|
| HRS | `HRS_USA/randhrs1992_2022v1.dta` | RAND HRS Longitudinal 2022 V1, May 2025 | `hhidpn` |
| HRS additional deficits | `HRS_USA/H_HRS_d.dta` | Harmonized HRS Ver.D | `hhidpn` |
| ELSA | `ELSA_UK/UKDA-5050-stata/stata/stata13_se/gh_elsa_h.dta` | Gateway Harmonized ELSA Ver.H, September 2025 | `idauniq` |
| SHARE | `SHARE_Europe/GH_SHARE_g.dta` | Gateway Harmonized SHARE Ver.G, SHARE release 9.0.0, October 2024 | `mergeid` |
| CHARLS | `CHARLS_China/H_CHARLS_D_Data.dta` | Harmonized CHARLS Ver.D, June 2021 | `ID` |
| MHAS | `MHAS_Mexico/H_MHAS_d.dta` | Gateway Harmonized MHAS Ver.D, May 2025 | `unhhidnp` |

Researchers obtain these released data through the cohorts' access processes. The [Gateway data-download guide](https://g2aging.org/additional-resources/gateway-blog/bbb35bc1-88a5-4826-b145-ce3e89b478fd) provides cohort routes. Harmonized SHARE additionally requires the licensed SHARE raw files and the corresponding Gateway creation code. The public repository includes the analysis code and source-variable mappings.

`mappings/baseline_variable_metadata.csv` records the actual available baseline fields and source paths relative to the configured data root. `mappings/complete_FI_encoding.csv` provides the exact retained fields. `mappings/incident_ADL_definition.csv` records the actual ADL aggregate variables. The data root is configurable; if the source file names differ, update the relative source mapping and source locations in `01_recover_baseline.R` consistently.

## Private intermediate objects

The numbered scripts create these objects in `MS5_RUNTIME_DIR`. The public package contains their aggregate results and schema descriptions.

1. `eligible_baseline.rds`: one row per eligible baseline participant, with unique `uid`, source ID, cohort, fixed baseline/follow-up wave numbers, age, clinical/context indicators, marital status, education, baseline `adl`, and `fu_adl`. Context is eight binary/missing indicators. HRS ID join preserves numeric identifier precision; IDs are formatted as non-scientific character strings.
2. `analysis/prepared_data.rds`: list with `main` and `panels`. `main` contains baseline/follow-up complete FI, item completeness, interview dates/status, follow-up interval, and inherited baseline/outcome fields. Each panel contains baseline/follow-up item-level 0–1 matrices, a retained-item map, and interview-date data. Panel rows match the corresponding cohort rows in `main`.
3. `analysis/results/candidate_question_bank.rds`: list containing `data`, `bank`, and `cohorts`. `data` retains baseline complete-FI participants and the 16 `q_` deficit columns.
4. `analysis/results/FI5_SHAP_holdout_<cohort>.rds`: the training-only latent-class model, class-specific and universal rankings, and training participant identifiers. The historical intermediate filename retains the original training workflow naming; the main reduction analysis derives scores for every length from the stored full rankings.
5. `analysis/results/FI_reduction_locked_scores.rds`: five held-out cohort lists with `uid`, reference FI, outcomes, class assignment, universal/phenotype FI1–FI16 matrices, complete-item matrices, locked training thresholds, and rankings.

Intermediate participant identifiers, fitted-model records containing training IDs, and item-level response matrices stay in the private runtime directory. Publication outputs are aggregate CSVs and figure files.

## Analysis sequence and verification scope

The README provides the sequential commands. Stage 01 recovers the baseline/follow-up backbone. Stage 02 constructs the complete FI using source metadata. Stage 03 constructs the 16-question bank. Stages 04–05 supply training latent classes and full SHAP rankings. Stage 06 creates locked score curves and training thresholds. Stages 07–10 evaluate performance, exploratory length rules, paired discrimination, and overlap sensitivity.

The original stored aggregate outputs accompany these scripts. The portable preparation audit evaluates stages 01–02 in memory and compares recovered participant membership and numeric fields with the frozen analysis objects. The audit intercepts all participant-level saves. Syntax and six-item boundary checks are separate. A licensed researcher can run the full chain using the listed releases and a private runtime directory; full model refitting is distinct from the preparation verification reported here.

## Remaining author-specific items

The author-supplied affiliations, degrees, contributions, funding status, and conflict status are recorded. Final submission still requires the corresponding author's email and the institution/documentation supporting the stated ethics exemption. The public repository address has been supplied by the authors; the root revision process handles its publication status and the final reproducibility statement.
