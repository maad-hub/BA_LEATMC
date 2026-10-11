# LEAMC analysis code

Analysis code of the bachelor thesis *Light Exposure Across the Menstrual Cycle* (Technical University of Munich).
The raw data are not included for data-protection reasons.

## Requirements

R 4.4.0 with the packages used for the thesis (versions of the final run):
LightLogR 0.10.3, menstrualcycleR 0.1.5 (`remotes::install_github("eisenlohrmoullab/menstrualcycleR")`), GGIR 3.3-8,
tidyverse 2.0.0, lme4 1.1-36, lmerTest 3.1-3, mgcv 1.9-4, nlme 3.1-169, psych 2.6.5, suncalc 0.5.3, rmcorr 0.7.0,
MASS 7.3-60.2, janitor, patchwork, corrplot, gt, ggbeeswarm.

## Running

1. Set the paths in `00_setup.R` (`BASE_DIR` = raw data, `OUTPUT_DIR` = results, `SYNC_OUTPUT_DIR`) and `GGIR_OUT`
   in `LEAMC_GGIR_sleep.R`.
2. Set the working directory to this folder and run `source("main.R")`.
   GGIR can take hours; with `RUN_GGIR <- FALSE` in `00_setup.R` the saved GGIR results are used.

## Scripts (in the order of `main.R`)

| Script | Content |
|---|---|
| `00_setup.R`, `00b_utils_io.R` | packages, paths, settings, helpers |
| `01_import_actlumus.R` | light-logger data (ActLumus) |
| `LEAMC_GGIR_sleep.R`, `02_import_actigraph.R`, `02b_import_activity_epochs.R` | actigraphy: sleep (GGIR) and activity time series |
| `03_import_redcap.R` | questionnaires (REDCap/MyCap), scoring, menses corrections |
| `04_join_trim.R`, `05_light_metrics.R`, `05b_activity_metrics.R` | daily data set, light and rest-activity metrics |
| `06_pacts.R` | cycle time (PACTS), phase labels, exclusions |
| `06b_phase_assignment.R` | extended phase assignment used for the main analyses |
| `07_descriptives.R` … `10b_thesis_numbers.R` | descriptives, missing data, LMMs and GAMMs, robustness checks, figures of the initial analysis, collected numbers |
| `revision.R` | Table 2, fixed-knot GAMMs, t-based CIs, study time, further robustness checks, rest-activity, sleep midpoint, repeated-measures correlations |
| `12_gamm_ar1.R` | final models for both phase assignments: GAMMs with AR(1), main tests P1–P4, Table 3, Table C1 (symptom index), adjusted GAMMs |
| `13_thesis_figures.R` | Figures 1–7 |
| `11_export.R` | export of data sets and workspace |
| `R/` | functions for 06b, 12 and 13 (`pa_` phase assignment, `ac_` model helpers, `e2_` final models and figures) |

In the code, `"A"` denotes the extended phase assignment (main analysis) and `"original"` the assignment based on
closed cycles only (sensitivity analysis).

## Where the results of the thesis come from (in `OUTPUT_DIR`)

| Thesis | File |
|---|---|
| Sample, Table 1, Cronbach's α, data per source, missing data, PACTS bins | `thesis_numbers.csv` |
| Phase assignment, study day by phase (both assignments) | `final_models/A_confound.csv`, `final_models/original_confound.csv`, `final_models/pacts_bins_by_assignment.csv` |
| Table 2 | `revision_output/table2_secondary_metrics.csv` (LMM, Holm over 16); GAMM with AR(1): `final_models/A_gamm.csv` (family `secondary`) |
| Table 3, Predictions 1–4 | `final_models/table3_both_assignments.csv`, `final_models/A_confirmatory.csv`, `final_models/A_lmm.csv` |
| Tables 4 and 5 (items) | `final_models/A_gamm.csv` (families `ls_items`, `drsp_items`, `impairment`; Holm in `p_holm_family`) |
| Table C1 | `final_models/A_robustness_csi.csv`, `final_models/adjusted_gamm.csv`, `robustness_table.csv`, `revision_output/results_for_placeholders.csv` (R15, R16, R20, R33, R35), `revision_output/R35_sensitivity_age_bmi.csv` |
| Table C2 (GAMMs without AR(1)) | `final_models/A_gamm.csv` (model `noAR`) |
| Other numbers in the text (IDs M*, R*) | `revision_output/results_for_placeholders.csv`, `revision_output/R*.csv` |
| Figures 1–7 | `final_models/figures/` |
| Figure B1 | `figures/fig_corr.png` |

## Analyses added after the initial analysis

As stated in the Methods: the extended phase assignment (`06b_phase_assignment.R`), GAMMs with an AR(1) residual
structure in which missing days remain gaps (`12_gamm_ar1.R`), and the GAMM as main test for Prediction 2.
The results of the initial analysis with the original assignment are part of the sensitivity analyses (`12_gamm_ar1.R`, variant `"original"`).
