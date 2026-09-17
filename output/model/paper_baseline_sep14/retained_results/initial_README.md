# Corrected initial equilibrium

Two exact numerical repetitions passed in job17655042. Loss: 179.2984242480. The nine structural search coordinates and empirical contract are unchanged from the selected point; the change is consistent float64 normalization of stored tenure probabilities in both Markov distribution operators. This is a verified candidate, not a converged optimizer.

Both repetitions have stationary nesting L1=2.2747454340e-14 and one-step nesting L1=2.6512812491e-13. Pension relative gap=3.0689413216e-12; property-tax rebate relative gap=3.8384144875e-7. The original gates remain unchanged. The two serialized checkpoint hashes differ, but the numerical signatures are exactly equal.

The largest weaknesses remain mean rooms and ownership at ages30–55. The correction changes the loss by less than0.000001; it does not establish an improved economic fit. Annual beta is at its enforced0.99 cap. Raw scorer bounds are retained in parameters_raw.csv.

## Complete target fit

| Moment | Target | Model | Gap | Weight | Loss contribution |
|---|---:|---:|---:|---:|---:|
| Initial model completed fertility | 2.1 | 2.100001 | 8.225487e-07 | — | — |
| Childless women, ages 40–44 | 0.1982788 | 0.1962584 | -0.002020369 | 35532.3 | 0.145039 |
| Exactly one child among mothers, ages 40–44 | 0.2136553 | 0.216439 | 0.002783641 | 26952.82 | 0.2088482 |
| Period mean first-birth age | 25.97626 | 26.09525 | 0.1189819 | 139.8281 | 1.979502 |
| First births at age 30+ | 0.249278 | 0.2321552 | -0.01712279 | 13866.07 | 4.06539 |
| Wealth / annual gross labor earnings | 6.145861 | 4.893427 | -1.252434 | 7.595098 | 11.91361 |
| Annual bequests / aggregate wealth | 0.0088 | 0.008582409 | -0.0002175905 | 5165289 | 0.2445539 |
| Old wealth/income p90 / median, ages 76–84 | 3.515935 | 4.4998 | 0.9838653 | 10.61636 | 10.27654 |
| Mean occupied rooms, capped at 9 | 5.561097 | 6.424652 | 0.8635551 | 128.0207 | 95.46854 |
| Ownership, heads 30–55 | 0.648334 | 0.5398962 | -0.1084379 | 2339.362 | 27.50802 |
| First-birth room response, −1 to +3 | 0.7202463 | 0.980806 | 0.2605598 | 137.5653 | 9.339497 |
| Rooms: 3+ versus 1–2 resident children (model dependent proxy) | 0.3470669 | 0.1720428 | -0.1750242 | 280.5281 | 8.593544 |
| Recent-parent ownership gap | 0.1628955 | 0.1441027 | -0.01879285 | 27055.82 | 9.555335 |

The first row is a separate normalization; twelve rows enter the scored objective. Full definitions, builders, samples, data vintages, uncertainty and weight provenance remain in target_fit.csv.

## Complete parameter and restriction table

| Parameter | Value | Lower | Upper | Near bound | Structural search coordinate |
|---|---:|---:|---:|---|---|
| beta_annual | 0.99 | 0.94 | 0.99 | True | True |
| kappa_fert | 0.3377342 | 0.02 | 50 | True | True |
| kappa_fert_continuation | 0.3978565 | 0.02 | 50 | True | True |
| chi | 1.049653 | 0.1 | 5 | False | True |
| H0 | 8.1121 | 0.2 | 80 | False | True |
| theta0 | 0.08105103 | 0 | 8 | False | True |
| theta1 | 0.08520357 | 0.02 | 16 | True | True |
| first_birth_fixed_cost | 0.2657648 | 0 | 8 | False | True |
| h_P | 2.3 | 0.1 | 2.3 | True | True |
| hbar_child_rooms | 0 | — | — | False | False |
| psi_child | 0.1489153 | — | — | False | False |
| payroll_tax | 0.179 | — | — | False | False |
| pension_period | 2.046361 | — | — | False | False |
| housing_supply_elasticity | 0.63 | — | — | False | False |
| tenure_choice_kappa | 0.005 | — | — | False | False |
| alpha_cons | 0.733 | — | — | False | False |
| sigma | 2 | — | — | False | False |

Bounds and near-bound flags follow the retained reporting rule, with the enforced annual-beta cap overlaid explicitly. Remaining nonstructural rows are normalizations, fiscal outcomes or externally fixed quantities; their full statuses and interpretations are in parameters.csv.

Source: corrected_initial_source_v2; frozen solver SHA256 2992412586b81cef3a3e58d92191bb51f54d3f9cc600d7675bbadaed7d1682da. Corrected full objective fingerprint4440ea07f4de957740ca6c04961d2806d9b9ef782c7a0e7dad4ce73e1db651b1. Only source-provenance fields and the numerical source-version label changed in the objective; the empirical-field fingerprint is unchanged.

## Native graph measurement check

The native housing graph reports uncapped demand 6.6398676277. The calibration
reports rooms capped at nine, 6.4246524480. The saved ten-room owner demand is
2.1521517966, corresponding to 0.2152151797 such owner households per total
household. Subtracting one room for each reconciles the two means to 1e-15;
see `rooms_cap_reconciliation.json`. This difference is not evidence of a
pre-choice/post-choice distribution mismatch. The actual scored solve uses
current realized tenure. Do not infer an income gradient from repeated plot
colors without checking the underlying age/type cells.
