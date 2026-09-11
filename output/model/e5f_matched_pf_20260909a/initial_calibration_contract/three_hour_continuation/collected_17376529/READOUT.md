# Completed initial-calibration refinement

**Job17376529 completed normally in54m35s.** All62 candidate jobs passed; two selected repetitions reproduce every numerical target row, every parameter row and loss exactly. The only target-row metadata difference is the checkpoint-file hash. All20 collected selected table/summary/PNG hashes checked.

The working objective falls from1180.224545 before both teaching-block jobs to869.413984 after the first batch and278.708079 after refinement (76.39% below the original point). The initial model fertility normalization remains2.1. This is pre-2007 initial calibration under the new utility and balanced pensions, not a2023 stationary recalibration, fitted historical path or policy benchmark. Same12scored targets plus one normalization; no changes to target definitions, weights, parameter bounds or numerical gates.

The three-hour limit was a ceiling. The fixed two-round search completed first; no further rounds were automatically launched.

## Complete target assessment

Gaps below are new model minus target, in each row’s native units. Weights include documented synthetic scales; the loss is a working minimum-distance criterion, not an optimal GMM statistic.

| Moment | Target | Before | New | New gap | Weight | Before loss | New loss |
|---|---:|---:|---:|---:|---:|---:|---:|
| Initial model completed fertility | 2.1 | 2.10003 | 2.10009 | 9.16457e-05 | — | — | — |
| Childless women, ages 40–44 | 0.198279 | 0.190789 | 0.217043 | 0.0187644 | 35532.3 | 1.99304 | 12.511 |
| Exactly one child among mothers, ages 40–44 | 0.213655 | 0.220384 | 0.190284 | -0.0233718 | 26952.8 | 1.22042 | 14.7227 |
| Period mean first-birth age | 25.9763 | 26.3139 | 26.002 | 0.0256865 | 139.828 | 15.9436 | 0.0922579 |
| First births at age 30+ | 0.249278 | 0.248208 | 0.225985 | -0.0232935 | 13866.1 | 0.0158684 | 7.52352 |
| Wealth / annual gross labor earnings | 6.14586 | 5.98092 | 5.92056 | -0.225302 | 7.5951 | 0.206629 | 0.385535 |
| Annual bequests / aggregate wealth | 0.0088 | 0.0114526 | 0.00851977 | -0.000280233 | 5.16529e+06 | 36.3431 | 0.405631 |
| Old wealth/income p90 / median, ages 76–84 | 3.51594 | 3.95692 | 4.34525 | 0.829319 | 10.6164 | 2.06456 | 7.30161 |
| Mean occupied rooms, capped at 9 | 5.5611 | 5.92039 | 6.22643 | 0.665336 | 128.021 | 16.5268 | 56.6712 |
| Ownership, heads 30–55 | 0.648334 | 0.604807 | 0.572675 | -0.0756586 | 2339.36 | 4.43222 | 13.391 |
| First-birth room response, −1 to +3 | 0.720246 | 0.528561 | 1.2659 | 0.545649 | 137.565 | 5.05459 | 40.9578 |
| Rooms: 3+ versus 1–2 resident children (model dependent proxy) | 0.347067 | 0.1864 | 0.166314 | -0.180753 | 280.528 | 7.24149 | 9.16535 |
| Recent-parent ownership gap | 0.162896 | -0.0377456 | 0.0975355 | -0.06536 | 27055.8 | 1089.18 | 115.58 |

## Complete parameters and restrictions

All nine structural search coordinates are present; psi is separately normalized. Fixed and derived parameters are retained below.

| Parameter | Estimate | Lower | Upper | Near bound? | Status |
|---|---:|---:|---:|---|---|
| beta_annual | 0.998944 | 0.94 | 0.9995 | True | diagnostic candidate; not a certified estimate |
| kappa_fert | 1.03226 | 0.02 | 50 | False | diagnostic candidate; not a certified estimate |
| kappa_fert_continuation | 1.01239 | 0.02 | 50 | False | diagnostic candidate; not a certified estimate |
| chi | 1.06631 | 0.1 | 5 | False | diagnostic candidate; not a certified estimate |
| H0 | 7.8693 | 0.2 | 80 | False | diagnostic candidate; not a certified estimate |
| theta0 | 0.189127 | 0 | 8 | False | diagnostic candidate; not a certified estimate |
| theta1 | 0.464858 | 0.02 | 16 | False | diagnostic candidate; not a certified estimate |
| first_birth_fixed_cost | 2.53751 | 0 | 8 | False | diagnostic candidate; not a certified estimate |
| h_P | 2.3 | 0.1 | 2.3 | True | diagnostic candidate; not a certified estimate |
| hbar_child_rooms | 0 | — | — | False | zero restriction |
| psi_child | 0.283219 | — | — | False | normalized to 2.1 |
| payroll_tax | 0.179 | — | — | False | externally fixed |
| pension_period | 2.04636 | — | — | False | budget derived |
| housing_supply_elasticity | 0.63 | — | — | False | externally fixed |
| tenure_choice_kappa | 0.005 | — | — | False | externally fixed |
| alpha_cons | 0.733 | — | — | False | externally fixed |
| sigma | 2 | — | — | False | externally fixed |

## Original diagnostic plots

All17 original selected plots are saved unchanged:

- [fertility_by_age](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/three_hour_continuation/collected_17376529/selected_standard_diagnostics/fertility_by_age.png)
- [fertility_policy_by_age_income_state](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/three_hour_continuation/collected_17376529/selected_standard_diagnostics/fertility_policy_by_age_income_state.png)
- [housing_by_age_income_state](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/three_hour_continuation/collected_17376529/selected_standard_diagnostics/housing_by_age_income_state.png)
- [housing_market](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/three_hour_continuation/collected_17376529/selected_standard_diagnostics/housing_market.png)
- [housing_prices](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/three_hour_continuation/collected_17376529/selected_standard_diagnostics/housing_prices.png)
- [income_state_outcomes](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/three_hour_continuation/collected_17376529/selected_standard_diagnostics/income_state_outcomes.png)
- [liquid_wealth_by_age_income_state](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/three_hour_continuation/collected_17376529/selected_standard_diagnostics/liquid_wealth_by_age_income_state.png)
- [market_clearing_by_market](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/three_hour_continuation/collected_17376529/selected_standard_diagnostics/market_clearing_by_market.png)
- [market_clearing_residuals](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/three_hour_continuation/collected_17376529/selected_standard_diagnostics/market_clearing_residuals.png)
- [owner_rungs](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/three_hour_continuation/collected_17376529/selected_standard_diagnostics/owner_rungs.png)
- [ownership_by_age](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/three_hour_continuation/collected_17376529/selected_standard_diagnostics/ownership_by_age.png)
- [ownership_by_age_income_state](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/three_hour_continuation/collected_17376529/selected_standard_diagnostics/ownership_by_age_income_state.png)
- [policy_childless_renter_age30](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/three_hour_continuation/collected_17376529/selected_standard_diagnostics/policy_childless_renter_age30.png)
- [policy_childless_renter_age42](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/three_hour_continuation/collected_17376529/selected_standard_diagnostics/policy_childless_renter_age42.png)
- [tenure_services](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/three_hour_continuation/collected_17376529/selected_standard_diagnostics/tenure_services.png)
- [wealth_dist_childless_renter_age30](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/three_hour_continuation/collected_17376529/selected_standard_diagnostics/wealth_dist_childless_renter_age30.png)
- [wealth_dist_childless_renter_age42](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/three_hour_continuation/collected_17376529/selected_standard_diagnostics/wealth_dist_childless_renter_age42.png)

## Remaining limitations

The recent-parent ownership observation remains the declared model approximation to the ACS definition. The parent ownership gap is9.75pp against16.29pp; first-birth housing response is1.266rooms against0.720rooms. Average rooms and overall ownership also remain imperfectly fitted. Exact reproduction and solver-gate passes do not establish a good overall fit, global identification, policy-function grid convergence or a certified historical/policy equilibrium.
