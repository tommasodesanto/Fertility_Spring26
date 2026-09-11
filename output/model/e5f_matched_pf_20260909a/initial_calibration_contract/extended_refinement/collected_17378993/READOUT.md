# Extended calibration: complete current readout

**Job17378993 completed normally in1h36m44s.** Allsix planned refinement rounds finished:178case evaluations,176verified and2housing-equilibrium rejections. The search continued past the rejections. The selected candidate was reproduced twice: all numerical target rows, parameter rows and loss match exactly.

Selected search candidate: `r5_joint_09`. Same working loss: **158.076391**, versus278.708079 at the start of this run and1180.224545 before the teaching-block searches. Improvement43.28% this run,86.61% across the sequence. All21collected selected artifact hashes verified; full objective independently recomputed.

**Scope:** pre-2007 stationary calibration with new utility and balanced pensions; twelve scored restrictions plus separate fertility normalization2.1. Same nine structural coordinates, targets, weights, bounds and solver gates. This is not a2023 stationary recalibration, fitted historical transition or policy benchmark. The older PDF contains an earlier candidate; the current full assessment is below.

**Economic assessment:** fertility levels fit reasonably; the recent-parent ownership gap is much closer to data. Overall ownership worsens and housing remains too large. The lower loss does not mean every moment improves. The housing-floor parameter has moved inside its upper bound; reported search-bound proximity flags remain for the two fertility shock scales andtheta1. These flags use the saved search convention and are not proof of binding constraints.

## Every target, weight and loss contribution

| Moment | Target | Previous | Current | Current gap | Weight | Previous loss | Current loss |
|---|---:|---:|---:|---:|---:|---:|---:|
| Initial model completed fertility | 2.1 | 2.10009 | 2.10015 | 0.000150605 | — | — | — |
| Childless women, ages 40–44 | 0.198279 | 0.217043 | 0.205018 | 0.00673908 | 35532.3 | 12.511 | 1.61371 |
| Exactly one child among mothers, ages 40–44 | 0.213655 | 0.190284 | 0.205322 | -0.00833344 | 26952.8 | 14.7227 | 1.87177 |
| Period mean first-birth age | 25.9763 | 26.002 | 26.1132 | 0.136904 | 139.828 | 0.0922579 | 2.62077 |
| First births at age 30+ | 0.249278 | 0.225985 | 0.230559 | -0.0187192 | 13866.1 | 7.52352 | 4.85881 |
| Wealth / annual gross labor earnings | 6.14586 | 5.92056 | 5.47059 | -0.675267 | 7.5951 | 0.385535 | 3.46326 |
| Annual bequests / aggregate wealth | 0.0088 | 0.00851977 | 0.00831932 | -0.000480683 | 5.16529e+06 | 0.405631 | 1.19347 |
| Old wealth/income p90 / median, ages 76–84 | 3.51594 | 4.34525 | 4.33505 | 0.819117 | 10.6164 | 7.30161 | 7.12307 |
| Mean occupied rooms, capped at 9 | 5.5611 | 6.22643 | 6.32135 | 0.760257 | 128.021 | 56.6712 | 73.9948 |
| Ownership, heads 30–55 | 0.648334 | 0.572675 | 0.534865 | -0.113469 | 2339.36 | 13.391 | 30.12 |
| First-birth room response, −1 to +3 | 0.720246 | 1.2659 | 1.11353 | 0.393287 | 137.565 | 40.9578 | 21.2779 |
| Rooms: 3+ versus 1–2 resident children (model dependent proxy) | 0.347067 | 0.166314 | 0.213083 | -0.133984 | 280.528 | 9.16535 | 5.03596 |
| Recent-parent ownership gap | 0.162896 | 0.0975355 | 0.149434 | -0.0134616 | 27055.8 | 115.58 | 4.90292 |

## Every parameter and restriction

| Parameter | Estimate | Lower | Upper | Near bound? | Status |
|---|---:|---:|---:|---|---|
| beta_annual | 0.995796 | 0.94 | 0.9995 | False | diagnostic candidate; not a certified estimate |
| kappa_fert | 0.334634 | 0.02 | 50 | True | diagnostic candidate; not a certified estimate |
| kappa_fert_continuation | 0.403731 | 0.02 | 50 | True | diagnostic candidate; not a certified estimate |
| chi | 1.04563 | 0.1 | 5 | False | diagnostic candidate; not a certified estimate |
| H0 | 8.2255 | 0.2 | 80 | False | diagnostic candidate; not a certified estimate |
| theta0 | 0.0895752 | 0 | 8 | False | diagnostic candidate; not a certified estimate |
| theta1 | 0.0770954 | 0.02 | 16 | True | diagnostic candidate; not a certified estimate |
| first_birth_fixed_cost | 0.293716 | 0 | 8 | False | diagnostic candidate; not a certified estimate |
| h_P | 2.1838 | 0.1 | 2.3 | False | diagnostic candidate; not a certified estimate |
| hbar_child_rooms | 0 | — | — | False | zero restriction |
| psi_child | 0.16526 | — | — | False | normalized to 2.1 |
| payroll_tax | 0.179 | — | — | False | externally fixed |
| pension_period | 2.04636 | — | — | False | budget derived |
| housing_supply_elasticity | 0.63 | — | — | False | externally fixed |
| tenure_choice_kappa | 0.005 | — | — | False | externally fixed |
| alpha_cons | 0.733 | — | — | False | externally fixed |
| sigma | 2 | — | — | False | externally fixed |

## Original diagnostic plots

All17selected original diagnostic plots are available below. Numerical reproducibility does not certify global policy accuracy or grid convergence. These show the new initial stationary candidate; they are not2023 transition policies.

- [fertility_by_age](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/extended_refinement/collected_17378993/selected_standard_diagnostics/fertility_by_age.png)
- [fertility_policy_by_age_income_state](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/extended_refinement/collected_17378993/selected_standard_diagnostics/fertility_policy_by_age_income_state.png)
- [housing_by_age_income_state](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/extended_refinement/collected_17378993/selected_standard_diagnostics/housing_by_age_income_state.png)
- [housing_market](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/extended_refinement/collected_17378993/selected_standard_diagnostics/housing_market.png)
- [housing_prices](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/extended_refinement/collected_17378993/selected_standard_diagnostics/housing_prices.png)
- [income_state_outcomes](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/extended_refinement/collected_17378993/selected_standard_diagnostics/income_state_outcomes.png)
- [liquid_wealth_by_age_income_state](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/extended_refinement/collected_17378993/selected_standard_diagnostics/liquid_wealth_by_age_income_state.png)
- [market_clearing_by_market](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/extended_refinement/collected_17378993/selected_standard_diagnostics/market_clearing_by_market.png)
- [market_clearing_residuals](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/extended_refinement/collected_17378993/selected_standard_diagnostics/market_clearing_residuals.png)
- [owner_rungs](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/extended_refinement/collected_17378993/selected_standard_diagnostics/owner_rungs.png)
- [ownership_by_age](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/extended_refinement/collected_17378993/selected_standard_diagnostics/ownership_by_age.png)
- [ownership_by_age_income_state](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/extended_refinement/collected_17378993/selected_standard_diagnostics/ownership_by_age_income_state.png)
- [policy_childless_renter_age30](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/extended_refinement/collected_17378993/selected_standard_diagnostics/policy_childless_renter_age30.png)
- [policy_childless_renter_age42](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/extended_refinement/collected_17378993/selected_standard_diagnostics/policy_childless_renter_age42.png)
- [tenure_services](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/extended_refinement/collected_17378993/selected_standard_diagnostics/tenure_services.png)
- [wealth_dist_childless_renter_age30](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/extended_refinement/collected_17378993/selected_standard_diagnostics/wealth_dist_childless_renter_age30.png)
- [wealth_dist_childless_renter_age42](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/extended_refinement/collected_17378993/selected_standard_diagnostics/wealth_dist_childless_renter_age42.png)

## Interpretation limits

The recent-parent model observer retains the documented approximation to the ACS household definition. Some weights use synthetic scales; this is a working minimum-distance criterion, not an optimal GMM statistic or a confidence test. No empirical target has been removed or reweighted. The initial fertility level2.1 is a separate author normalization rather than an empirical period TFR fit. Historical preference-path fitting, horizon certification, updated policy comparisons and the separate matched2023 static recalibration remain outstanding. No automatic next run has been launched.
