# Final income search 18049121 readout

Mechanical readout of the terminal verified selection: 96 proposals, 89 valid and 7 rejected; selected case 60 has objective/loss `353.6588729140903`. The selected point has two exact native repetitions. The search seed `psi=0.23949950404168222` is distinct from the normalized parameter `psi_child=0.19227029432243803` in the 17-parameter table.

## Target fit

| Target | Target value | Model | Gap | Actual weight | Loss contribution |
|---|---:|---:|---:|---:|---:|
| Initial model completed fertility | 2.1 | 2.10027 | 0.000270543 | — | — |
| Childless women, ages 40–44 | 0.198279 | 0.225169 | 0.0268907 | 35532.3 | 25.6938 |
| Exactly one child among mothers, ages 40–44 | 0.213655 | 0.197956 | -0.0156996 | 26952.8 | 6.64326 |
| Period mean first-birth age | 25.9763 | 26.7001 | 0.72383 | 139.828 | 73.26 |
| First births at age 30+ | 0.249278 | 0.259829 | 0.0105508 | 13866.1 | 1.54357 |
| Wealth / annual gross labor earnings | 6.14586 | 6.72243 | 0.576571 | 7.5951 | 2.52487 |
| Annual bequests / aggregate wealth | 0.0088 | 0.0079772 | -0.000822802 | 5.16529e+06 | 3.49691 |
| Old wealth/income p90 / median, ages 76–84 | 3.51594 | 4.28269 | 0.766753 | 10.6164 | 6.24147 |
| Mean occupied rooms, capped at 9 | 5.5611 | 6.67698 | 1.11588 | 128.021 | 159.411 |
| Ownership, heads 30–55 | 0.648334 | 0.484271 | -0.164063 | 2339.36 | 62.968 |
| First-birth room response, −1 to +3 | 0.720246 | 1.01391 | 0.293668 | 137.565 | 11.8638 |
| Rooms: 3+ versus 1–2 resident children (model dependent proxy) | 0.347067 | 0.344338 | -0.00272857 | 280.528 | 0.00208856 |
| Recent-parent ownership gap | 0.162896 | 0.162282 | -0.000613321 | 27055.8 | 0.0101774 |

The 12 scored loss contributions sum to `353.6588729140903`.

## Parameters

The nine structural active bounds below are taken from `overnight/plan.remote.json#parameter_bounds`, with the active beta upper bound set to `0.99`. Near-bound flags use the active bounds; the remaining rows report their restriction or external status.

| Parameter | Estimate | Active lower | Active upper | Restriction / external status | Near bound |
|---|---:|---:|---:|---|---|
| beta_annual | 0.99 | 0.94 | 0.99 | diagnostic candidate; not a certified estimate | true |
| kappa_fert | 0.197771 | 0.02 | 50 | diagnostic candidate; not a certified estimate | true |
| kappa_fert_continuation | 0.346494 | 0.02 | 50 | diagnostic candidate; not a certified estimate | true |
| chi | 0.953059 | 0.1 | 5 | diagnostic candidate; not a certified estimate | false |
| H0 | 10.1554 | 0.2 | 80 | diagnostic candidate; not a certified estimate | false |
| theta0 | 0.170783 | 0 | 8 | diagnostic candidate; not a certified estimate | false |
| theta1 | 0.0734872 | 0.02 | 16 | diagnostic candidate; not a certified estimate | true |
| first_birth_fixed_cost | 0.167382 | 0 | 8 | diagnostic candidate; not a certified estimate | false |
| h_P | 2.3 | 0.1 | 2.3 | diagnostic candidate; not a certified estimate | true |
| hbar_child_rooms | 0 | — | — | zero restriction | false |
| psi_child | 0.19227 | — | — | normalized to 2.1 | false |
| payroll_tax | 0.179 | — | — | externally fixed | false |
| pension_period | 2.04636 | — | — | budget derived | false |
| housing_supply_elasticity | 0.63 | — | — | externally fixed | false |
| tenure_choice_kappa | 0.005 | — | — | externally fixed | false |
| alpha_cons | 0.733 | — | — | externally fixed | false |
| sigma | 2 | — | — | externally fixed | false |

Verification receipt flags are `numeric_fit_equal=true`, `exact_objective=true`, and `native_exact_loss_equality=true`. The packet contains 13 fit rows, 17 parameter rows, and 17 selected diagnostic PNGs; their SHA-256 values are recorded in `receipt.json`.

Raw parameter CSV preserves the scorer’s generic beta upper bound 0.9995; the table above uses the actual search restriction 0.99. This is a finite local multivariate diagnostic, not a converged or adopted calibration.
