# Follow-up availability and sensitivity analysis

## Scope and selection stages

Analyses use the existing five-cohort, baseline complete-FI sample and the locked common-score risk sets. Selection is conditional on the inherited, same-baseline-wave eligibility criteria. Baseline ages range from 65 years upward in each cohort; the raw-data crosscheck matches the prepared selected-wave interview status for every participant.

| Cohort | Baseline complete FI | Common scores | Common-score exclusion | Follow-up FI unavailable within common | Retained |
|---|---:|---:|---:|---:|---:|
| HRS | 8,349 | 8,075 | 274 | 1,565 | 6,510 |
| ELSA | 4,197 | 4,196 | 1 | 1,615 | 2,581 |
| SHARE | 12,560 | 12,552 | 8 | 6,701 | 5,851 |
| CHARLS | 4,021 | 3,946 | 75 | 1,016 | 2,930 |
| MHAS | 5,827 | 5,433 | 394 | 1,444 | 3,989 |

## Retained and unretained comparison

Each baseline variable is summarized within its observed subset; separate missing counts accompany means or proportions. Continuous standardized mean differences use the square root of the mean of the retained and comparison-group variances. The sign is retained minus unretained. Seven-condition burden sums hypertension, diabetes, heart disease, cancer, stroke, arthritis and lung disease among participants with all seven conditions observed.

| SHARE variable | Retained | Unretained | SMD |
|---|---:|---:|---:|
| Age, years | 73.036 | 75.498 | -0.378 |
| Female | 0.573 | 0.543 | 0.060 |
| Baseline complete FI | 0.205 | 0.263 | -0.384 |
| Seven-condition burden | 1.849 | 1.936 | -0.089 |
| Married or partnered | 0.677 | 0.633 | 0.093 |
| Low education | 0.601 | 0.638 | -0.077 |
| Any baseline ADL limitation | 0.136 | 0.237 | -0.262 |
| Any baseline IADL limitation | 0.099 | 0.250 | -0.404 |

## Source-labelled mortality and other non-observation

The original Stata variable and value labels identify each selected-wave interview-status code. The status fields are complete in all five cohorts and match the prepared fields. Recorded death codes are HRS 2/3/5/6, ELSA 3/5/6, and SHARE/CHARLS/MHAS 5/6. Code 9 identifies unresolved alive/dead status; code 10 in SHARE identifies country absence from the wave. Labels containing both alive and died are classified as unresolved status.

Among the 6,701 SHARE common-score participants whose follow-up FI is unavailable, 1,577 have recorded-death status, 2,502 are explicitly alive nonrespondents, 1,161 have unresolved vital status, 1,334 belong to a country absent from the wave, and 127 are alive respondents with an insufficiently observed FI. The eight common-score exclusions are tabulated separately.

At baseline-complete-FI level, SHARE country-wave-absence codes comprise Greece (952), Ireland (246) and Israel (139). Israel additionally contains 895 explicitly alive nonrespondents and zero follow-up FI observations; 893 meet the common-score criterion.

ELSA has 2,651 alive/respondent statuses and 1,546 unresolved vital statuses. An ancillary local end-of-life linkage identifies 103 ELSA code-9 participants, with death year observed for 102; 52 deaths in 2018–2020 predate the 2021 start of wave 10. Another 50 have death years in 2021–2023, and one death year is unobserved. The linkage covers a subset of original unresolved cases and is reported as ancillary date evidence. SHARE ancillary end-of-life data identify three code-9 deaths in 2006–2009, predating its wave-4 start. Primary status groups retain their selected-wave source definitions.

## IPW methods paragraph

We examined outcome-availability sensitivity among participants with common locked scores, baseline complete FI, and an explicitly alive follow-up status. The target was further restricted to countries with at least one observed follow-up FI. Country strata with zero observed follow-up FI were reported separately. Cohort-specific logistic models estimated follow-up FI observability from baseline age, sex, complete FI, seven-condition burden, marital status, low education, baseline ADL and IADL limitation, walking difficulty, and self-rated health; age and FI used three-degree-of-freedom natural splines. SHARE models included source-derived country fixed effects. Missing baseline covariates received cohort-specific median imputation with missing indicators. Stabilized inverse-probability weights used the cohort observability fraction as numerator, fitted probabilities clipped to 0.05–0.99, and retained-person weights truncated at their first and 99th percentiles. We estimated weighted AUROC for follow-up FI ≥0.25 and weighted mid-distribution rank correlation with continuous follow-up FI. Percentile 95% intervals used 200 participant bootstraps with the observability model refitted in each replicate. Weight ranges, clipping counts, truncation counts, effective sample sizes, and baseline balance are supplied.

This sensitivity analysis addresses covariate-associated outcome availability within the explicitly alive, country-supported population. Recorded mortality, unresolved vital status and zero-observation country strata define separate estimand boundaries. ELSA’s verified-alive target primarily addresses FI-item observability among interviewed participants. The five cohorts use their original selected baseline participants.

## IPW results

| Cohort | Target / observed | FI6 AUROC: observed → IPW (95% CI) | Full FI AUROC: observed → IPW (95% CI) | Effective sample size |
|---|---:|---|---|---:|
| HRS | 7,317 / 6,510 | 0.868 → 0.870 (0.861–0.878) | 0.918 → 0.919 (0.913–0.926) | 6484 |
| ELSA | 2,650 / 2,581 | 0.853 → 0.855 (0.837–0.869) | 0.888 → 0.890 (0.876–0.903) | 2578 |
| SHARE | 7,587 / 5,851 | 0.810 → 0.812 (0.801–0.822) | 0.849 → 0.850 (0.840–0.860) | 5786 |
| CHARLS | 3,344 / 2,930 | 0.797 → 0.800 (0.784–0.815) | 0.842 → 0.843 (0.831–0.858) | 2886 |
| MHAS | 4,557 / 3,989 | 0.793 → 0.797 (0.781–0.811) | 0.844 → 0.846 (0.833–0.857) | 3964 |

## Shared five-task ADL methods paragraph

A task-harmonized sensitivity endpoint used bathing, dressing, eating, getting in or out of bed, and toileting in every cohort. Baseline freedom from limitation required all five task responses observed and all five scored zero. At follow-up, any observed limitation established a positive endpoint; a negative endpoint required all five observed responses scored zero. Analyses used the same locked common-score requirement and evaluated incident limitation among participants free of all five limitations at baseline. AUROC differences and above-chance AUROC retention, (AUROC_FI6 − 0.5)/(AUROC_full FI − 0.5), used paired resampling with 200 participant bootstraps. The main cohort-specific endpoint results remain the primary analysis.

| Cohort | Shared endpoint n / events | FI6 AUROC | Full FI AUROC | Above-chance AUROC retention (95% CI) |
|---|---:|---:|---:|---|
| HRS | 5,422 / 587 | 0.768 | 0.784 | 94.1% (90.3–98.1%) |
| ELSA | 2,181 / 355 | 0.738 | 0.755 | 93.5% (86.7–100.2%) |
| SHARE | 5,130 / 833 | 0.685 | 0.704 | 90.8% (84.5–96.7%) |
| CHARLS | 2,288 / 464 | 0.680 | 0.705 | 87.7% (78.5–96.4%) |
| MHAS | 3,467 / 763 | 0.667 | 0.701 | 83.1% (76.0–88.6%) |

The task-harmonized analysis contains 18,488 participants. Mean cohort-specific above-chance AUROC retention is 89.8%. Task-level source variables, coding rules, missing counts, sample flows, and primary/shared endpoint comparisons accompany the result table.

## Release and reproducibility

Public aggregate results are the CSV tables, this report, field-gap notes, and R scripts. Participant-level RDS audits and the data-structure inventory contain local analysis information and are for the controlled workspace. Original source files are read-only. The original data preparation and score files are required to reproduce the scripts; source licences determine participant-level sharing.

