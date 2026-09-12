# Completed capped-beta calibration — September 11, evening

Loss **159.238986075**, versus 158.076390747 for the unrestricted reference: 0.735% higher under identical targets and weights. Annual beta was estimated in [0.94,0.99] and the best tested point is at 0.99. All nine structural parameters were free. This is a provisional initial-state working minimum-distance fit, not a certified historical calibration or identified global optimum.

The bounded search completed 84 new trials and two exact final repetitions. One trial failed its strict mass gate and was excluded; all selected-point gates passed. Seventeen standard diagnostic graphs were produced. The separate fixed-beta comparison stopped after a recurring mass-gate failure; it is not the selected result.

## Full target fit

| Moment | Target | Model | Gap | Weight | Loss contribution |
|---|---:|---:|---:|---:|---:|
| Initial model completed fertility | 2.1 | 2.100001 | 1.024836e-06 | — | — |
| Childless women, ages 40–44 | 0.1982788 | 0.2062257 | 0.007946926 | 35532.3 | 2.243994 |
| Exactly one child among mothers, ages 40–44 | 0.2136553 | 0.2069008 | -0.006754562 | 26952.82 | 1.229698 |
| Period mean first-birth age | 25.97626 | 26.09785 | 0.1215872 | 139.8281 | 2.067139 |
| First births at age 30+ | 0.249278 | 0.2317489 | -0.0175291 | 13866.07 | 4.260615 |
| Wealth / annual gross labor earnings | 6.145861 | 4.847472 | -1.29839 | 7.595098 | 12.80394 |
| Annual bequests / aggregate wealth | 0.0088 | 0.008552889 | -0.0002471111 | 5165289 | 0.3154126 |
| Old wealth/income p90 / median, ages 76–84 | 3.515935 | 4.434562 | 0.918627 | 10.61636 | 8.958888 |
| Mean occupied rooms, capped at 9 | 5.561097 | 6.288719 | 0.7276215 | 128.0207 | 67.77839 |
| Ownership, heads 30–55 | 0.648334 | 0.5229541 | -0.12538 | 2339.362 | 36.77509 |
| First-birth room response, −1 to +3 | 0.7202463 | 1.013756 | 0.2935095 | 137.5653 | 11.85095 |
| Rooms: 3+ versus 1–2 resident children (model dependent proxy) | 0.3470669 | 0.1965613 | -0.1505057 | 280.5281 | 6.354509 |
| Recent-parent ownership gap | 0.1628955 | 0.1498559 | -0.01303963 | 27055.82 | 4.600358 |

The 2.1 completed-fertility value is a separately imposed model normalization, not an observed female period TFR. CPS age interpolation and the ACS recent-parent mapping retain their declared approximations. All row definitions and provenance are preserved in selected_capped_score.json; weights are the frozen working weights, not a claim of a fully certified empirical SMM specification.

## Full parameter and restriction table

| Parameter | Estimate | Lower | Upper | Near bound | Status |
|---|---:|---:|---:|---|---|
| beta_annual | 0.99 | 0.94 | 0.99 | True | estimated_capped |
| kappa_fert | 0.2971521 | 0.02 | 50 | True | diagnostic candidate; not a certified estimate |
| kappa_fert_continuation | 0.363372 | 0.02 | 50 | True | diagnostic candidate; not a certified estimate |
| chi | 1.047576 | 0.1 | 5 | False | diagnostic candidate; not a certified estimate |
| H0 | 8.185325 | 0.2 | 80 | False | diagnostic candidate; not a certified estimate |
| theta0 | 0.08354497 | 0 | 8 | False | diagnostic candidate; not a certified estimate |
| theta1 | 0.05167858 | 0.02 | 16 | True | diagnostic candidate; not a certified estimate |
| first_birth_fixed_cost | 0.154874 | 0 | 8 | False | diagnostic candidate; not a certified estimate |
| h_P | 2.3 | 0.1 | 2.3 | True | diagnostic candidate; not a certified estimate |
| hbar_child_rooms | 0 | — | — | False | zero restriction |
| psi_child | 0.1612272 | — | — | False | normalized to 2.1 |
| payroll_tax | 0.179 | — | — | False | externally fixed |
| pension_period | 2.046361 | — | — | False | budget derived |
| housing_supply_elasticity | 0.63 | — | — | False | externally fixed |
| tenure_choice_kappa | 0.005 | — | — | False | externally fixed |
| alpha_cons | 0.733 | — | — | False | externally fixed |
| sigma | 2 | — | — | False | externally fixed |

Beta and the parent housing floor h_P are exactly at their upper bounds. Other near-bound flags use the existing 1% of range convention and do not establish a binding optimum.

## Historical and policy status at the evening check

The following transition trials use the earlier unrestricted initial calibration, not the newly selected capped candidate. All three six-date roots passed market, fiscal, exact replay and checkpoint gates. No 28-date path passed its finite-equilibrium or terminal-horizon checks. No preference path has been adaptively estimated, and the household-based fertility diagnostic still needs an agreed mapping to female period TFR. Policy results under the new utility/pension/calibration combination are not verified.

| Trial preference decline | Maximum housing gap | Maximum pension gap | Finite equilibrium passes | Horizon passes |
|---|---:|---:|---|---|
| -0.025 | 0.80922% | 0.023208% | False | False |
| -0.05 | 1.10561% | 0.028509% | False | False |
| -0.1 | 1.48018% | 0.107444% | False | False |

Required residual tolerances are 0.02% for housing and 0.0001% for pensions. These are finite-path root failures after the budgeted iterations, not evidence that the budget identity was removed. The terminal distribution checks also fail. All cluster jobs are finished; no new run was submitted during this assessment.

## Decision proposed to the author

Provisionally freeze this initial candidate and the sequential/new-utility/balanced-pension specification. Concentrate further work on the empirical fertility observation mapping, the fitted announced preference path, and convergence/horizon checks before drawing policy conclusions. Do not use old-policy numbers as results for this revised calibration.

## Existing diagnostic plots

- [fertility_by_age](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/fixed_beta_profiles/collected_17425504/selected_standard_diagnostics/fertility_by_age.png)
- [fertility_policy_by_age_income_state](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/fixed_beta_profiles/collected_17425504/selected_standard_diagnostics/fertility_policy_by_age_income_state.png)
- [housing_by_age_income_state](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/fixed_beta_profiles/collected_17425504/selected_standard_diagnostics/housing_by_age_income_state.png)
- [housing_market](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/fixed_beta_profiles/collected_17425504/selected_standard_diagnostics/housing_market.png)
- [housing_prices](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/fixed_beta_profiles/collected_17425504/selected_standard_diagnostics/housing_prices.png)
- [income_state_outcomes](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/fixed_beta_profiles/collected_17425504/selected_standard_diagnostics/income_state_outcomes.png)
- [liquid_wealth_by_age_income_state](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/fixed_beta_profiles/collected_17425504/selected_standard_diagnostics/liquid_wealth_by_age_income_state.png)
- [market_clearing_by_market](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/fixed_beta_profiles/collected_17425504/selected_standard_diagnostics/market_clearing_by_market.png)
- [market_clearing_residuals](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/fixed_beta_profiles/collected_17425504/selected_standard_diagnostics/market_clearing_residuals.png)
- [owner_rungs](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/fixed_beta_profiles/collected_17425504/selected_standard_diagnostics/owner_rungs.png)
- [ownership_by_age](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/fixed_beta_profiles/collected_17425504/selected_standard_diagnostics/ownership_by_age.png)
- [ownership_by_age_income_state](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/fixed_beta_profiles/collected_17425504/selected_standard_diagnostics/ownership_by_age_income_state.png)
- [policy_childless_renter_age30](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/fixed_beta_profiles/collected_17425504/selected_standard_diagnostics/policy_childless_renter_age30.png)
- [policy_childless_renter_age42](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/fixed_beta_profiles/collected_17425504/selected_standard_diagnostics/policy_childless_renter_age42.png)
- [tenure_services](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/fixed_beta_profiles/collected_17425504/selected_standard_diagnostics/tenure_services.png)
- [wealth_dist_childless_renter_age30](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/fixed_beta_profiles/collected_17425504/selected_standard_diagnostics/wealth_dist_childless_renter_age30.png)
- [wealth_dist_childless_renter_age42](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/fixed_beta_profiles/collected_17425504/selected_standard_diagnostics/wealth_dist_childless_renter_age42.png)
