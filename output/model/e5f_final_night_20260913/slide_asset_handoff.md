# Slide asset handoff (provisional; pending 12:30 EDT freeze)

Inventory verified against `README.md`, the two figure manifests, validation CSV headers, and `verified_history_readout_qa.json`.

## Main report and page map

- Main 41-page report: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_final_night_20260913/verified_history_readout.pdf`
- Pages 2–8: complete initial target fit, historical fit, untargeted 2023 validation, paired tax tables, and parameter restrictions.
- Pages 9–35: all 17 native diagnostics for initial equilibrium, baseline policy, and higher-tax policy.
- Pages 36–41: five historical/2023 comparison figures and supplemental policy comparison.
- QA: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_final_night_20260913/verified_history_readout_qa.json` (verified; 41 pages).

## Historical/data comparison figures (five PDF/PNG pairs)

- Historical fertility: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_final_night_20260913/history_A0_6/figures/historical_fertility.pdf` and `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_final_night_20260913/history_A0_6/figures/historical_fertility.png`
- Prices and quantities path: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_final_night_20260913/history_A0_6/figures/prices_quantities_path.pdf` and `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_final_night_20260913/history_A0_6/figures/prices_quantities_path.png`
- 2023 equilibrium: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_final_night_20260913/history_A0_6/figures/equilibrium_2023.pdf` and `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_final_night_20260913/history_A0_6/figures/equilibrium_2023.png`
- 2023 lifecycle: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_final_night_20260913/history_A0_6/figures/lifecycle_2023.pdf` and `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_final_night_20260913/history_A0_6/figures/lifecycle_2023.png`
- 2023 fertility by age: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_final_night_20260913/history_A0_6/figures/fertility_age_2023.pdf` and `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_final_night_20260913/history_A0_6/figures/fertility_age_2023.png`

## Supplemental paired policy comparison

- `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_final_night_20260913/history_A0_6/policy_comparison/policy_comparison.pdf`
- `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_final_night_20260913/history_A0_6/policy_comparison/policy_comparison.png`

## Numerical source CSVs (headers verified)

- 13-row initial fit: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_final_night_20260913/corrected_initial/target_fit.csv` (header begins `actual_weight,...`).
- 17-parameter table: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_final_night_20260913/corrected_initial/parameters.csv` (header begins `estimate,...`).
- 13-row untargeted 2023 validation: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_final_night_20260913/history_A0_6/validation/validation_2023.csv` (header `moment_key,moment,data,model,gap,data_vintage,model_year,decimals,measurement_note`).
- 66-row policy comparison: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_final_night_20260913/history_A0_6/policy_comparison/comparison.csv` (header `decision_year,fertility_window_end,metric,baseline,reform,absolute_change,percent_change,percentage_point_change`).

## Established limitations (verbatim from README)

- “The six-period results remain **provisional**: horizon adequacy is unverified, and production eligibility remains false.”
- “Actual empirical vintages and measurement qualifications remain explicit.”
- “A0 removes all post-2023 migration; A+ retains the supplied migration sensitivity.”
- “These are dated flow/level effects, not cumulative births or long-run stationary comparisons.”
- “Both 1% and 2% annual property taxes return revenue equally per current head.”
- “Every accepted root passes the finite housing, PAYGO, equal-rebate and exact-replay gates.”

Handoff remains provisional pending the 12:30 EDT freeze.
