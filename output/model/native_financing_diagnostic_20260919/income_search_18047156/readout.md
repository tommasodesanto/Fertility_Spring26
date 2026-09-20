# Income search 18047156 receipt

Mechanical readout of the verified selection (case 3; objective/loss `502.74561411262744`). The source run completed 16 proposals and the selected point was repeated twice.

## Target fit

| Target | Target value | Model | Gap | Actual weight | Loss contribution |
|---|---:|---:|---:|---:|---:|
| Initial model completed fertility | 2.1 | 2.10018 | 0.000175118 | — | — |
| Childless women, ages 40–44 | 0.198279 | 0.22236 | 0.0240812 | 35532.3 | 20.6053 |
| Exactly one child among mothers, ages 40–44 | 0.213655 | 0.203773 | -0.00988243 | 26952.8 | 2.63228 |
| Period mean first-birth age | 25.9763 | 27.127 | 1.15074 | 139.828 | 185.161 |
| First births at age 30+ | 0.249278 | 0.285435 | 0.0361569 | 13866.1 | 18.1274 |
| Wealth / annual gross labor earnings | 6.14586 | 6.69245 | 0.546593 | 7.5951 | 2.26914 |
| Annual bequests / aggregate wealth | 0.0088 | 0.00701284 | -0.00178716 | 5.16529e+06 | 16.4976 |
| Old wealth/income p90 / median, ages 76–84 | 3.51594 | 4.54856 | 1.03262 | 10.6164 | 11.3203 |
| Mean occupied rooms, capped at 9 | 5.5611 | 6.18288 | 0.621779 | 128.021 | 49.494 |
| Ownership, heads 30–55 | 0.648334 | 0.489777 | -0.158557 | 2339.36 | 58.8123 |
| First-birth room response, −1 to +3 | 0.720246 | 1.30832 | 0.588074 | 137.565 | 47.5743 |
| Rooms: 3+ versus 1–2 resident children (model dependent proxy) | 0.347067 | 0.339694 | -0.00737341 | 280.528 | 0.0152515 |
| Recent-parent ownership gap | 0.162896 | 0.105144 | -0.0577512 | 27055.8 | 90.2367 |

## Parameters

Active bounds for the nine structural coordinates are copied from `earnings_candidate/search_plan.remote.json`; the active-bound flag uses the run’s near-bound row, with the explicit beta upper coordinate recorded as `0.99`.

| Parameter | Estimate | Active lower | Active upper | Restriction / external status | Near bound |
|---|---:|---:|---:|---|---|
| beta_annual | 0.99 | 0.94 | 0.99 | diagnostic candidate; not a certified estimate | true |
| kappa_fert | 0.337734 | 0.02 | 50 | diagnostic candidate; not a certified estimate | true |
| kappa_fert_continuation | 0.477428 | 0.02 | 50 | diagnostic candidate; not a certified estimate | true |
| chi | 1.04965 | 0.1 | 5 | diagnostic candidate; not a certified estimate | false |
| H0 | 8.1121 | 0.2 | 80 | diagnostic candidate; not a certified estimate | false |
| theta0 | 0.081051 | 0 | 8 | diagnostic candidate; not a certified estimate | false |
| theta1 | 0.0852036 | 0.02 | 16 | diagnostic candidate; not a certified estimate | true |
| first_birth_fixed_cost | 0.265765 | 0 | 8 | diagnostic candidate; not a certified estimate | false |
| h_P | 2.3 | 0.1 | 2.3 | diagnostic candidate; not a certified estimate | true |
| hbar_child_rooms | 0 | — | — | zero restriction | false |
| psi_child | 0.2395 | — | — | normalized to 2.1 | false |
| payroll_tax | 0.179 | — | — | externally fixed | false |
| pension_period | 2.04636 | — | — | budget derived | false |
| housing_supply_elasticity | 0.63 | — | — | externally fixed | false |
| tenure_choice_kappa | 0.005 | — | — | externally fixed | false |
| alpha_cons | 0.733 | — | — | externally fixed | false |
| sigma | 2 | — | — | externally fixed | false |

Verification receipt: `numeric_fit_equal=true`, `exact_objective=true`, `native_exact_loss_equality=true`; 13 fit rows, 17 parameter rows, and 17 native diagnostic PNGs were collected.
