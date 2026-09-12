# Completed overnight initial calibration

Loss **158.541191063**, versus 159.238986075 at the start: **0.438% lower**. Three bounded batches completed 252 new search trials. Each selected batch winner passed two exact repetitions. This remains a provisional initial calibration, not a fitted historical path or policy result.

Beta remains estimated and capped at 0.99. The complete target/weight fingerprint is unchanged. Initial fertility 2.1 is a separate normalization. All 17 standard graphs are in `selected_standard_diagnostics/`. The historical workers retain their previously pinned calibration.

| Moment | Target | Model | Gap | Weight | Loss contribution |
|---|---:|---:|---:|---:|---:|
| Initial model completed fertility | 2.1 | 2.1 | 8.9827e-07 | — | — |
| Childless women, ages 40–44 | 0.198279 | 0.206956 | 0.00867741 | 35532.3 | 2.67549 |
| Exactly one child among mothers, ages 40–44 | 0.213655 | 0.20658 | -0.00707496 | 26952.8 | 1.34912 |
| Period mean first-birth age | 25.9763 | 26.1032 | 0.126914 | 139.828 | 2.25224 |
| First births at age 30+ | 0.249278 | 0.232027 | -0.0172509 | 13866.1 | 4.12646 |
| Wealth / annual gross labor earnings | 6.14586 | 4.84892 | -1.29695 | 7.5951 | 12.7755 |
| Annual bequests / aggregate wealth | 0.0088 | 0.00858255 | -0.000217451 | 5.16529e+06 | 0.24424 |
| Old wealth/income p90 / median, ages 76–84 | 3.51594 | 4.43839 | 0.922459 | 10.6164 | 9.03378 |
| Mean occupied rooms, capped at 9 | 5.5611 | 6.29059 | 0.729492 | 128.021 | 68.1272 |
| Ownership, heads 30–55 | 0.648334 | 0.52444 | -0.123894 | 2339.36 | 35.9088 |
| First-birth room response, −1 to +3 | 0.720246 | 1.01103 | 0.29078 | 137.565 | 11.6315 |
| Rooms: 3+ versus 1–2 resident children (model dependent proxy) | 0.347067 | 0.200474 | -0.146593 | 280.528 | 6.02841 |
| Recent-parent ownership gap | 0.162896 | 0.15016 | -0.0127356 | 27055.8 | 4.38837 |

| Parameter/restriction | Estimate | Lower | Upper | Near bound | Status |
|---|---:|---:|---:|---|---|
| beta_annual | 0.99 | 0.94 | 0.99 | True | estimated_capped |
| kappa_fert | 0.279392 | 0.02 | 50 | True | diagnostic candidate; not a certified estimate |
| kappa_fert_continuation | 0.346867 | 0.02 | 50 | True | diagnostic candidate; not a certified estimate |
| chi | 1.04803 | 0.1 | 5 | False | diagnostic candidate; not a certified estimate |
| H0 | 8.19051 | 0.2 | 80 | False | diagnostic candidate; not a certified estimate |
| theta0 | 0.0869545 | 0 | 8 | False | diagnostic candidate; not a certified estimate |
| theta1 | 0.0571137 | 0.02 | 16 | True | diagnostic candidate; not a certified estimate |
| first_birth_fixed_cost | 0.1268 | 0 | 8 | False | diagnostic candidate; not a certified estimate |
| h_P | 2.3 | 0.1 | 2.3 | True | diagnostic candidate; not a certified estimate |
| hbar_child_rooms | 0 | — | — | False | zero restriction |
| psi_child | 0.159621 | — | — | False | normalized to 2.1 |
| payroll_tax | 0.179 | — | — | False | externally fixed |
| pension_period | 2.04636 | — | — | False | budget derived |
| housing_supply_elasticity | 0.63 | — | — | False | externally fixed |
| tenure_choice_kappa | 0.005 | — | — | False | externally fixed |
| alpha_cons | 0.733 | — | — | False | externally fixed |
| sigma | 2 | — | — | False | externally fixed |

The 12 scored contributions recompute the loss exactly. The full CSVs preserve empirical provenance and uncertainty. Derived fiscal entries in a parameter table are not dated transition pension balances; consult the actual fiscal receipts for solved-state balances.
