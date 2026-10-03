# Matched-wealth local calibration: ten-chain collection

Collected October 3, 2026. All **10/10** chains passed fresh native selected-point verification. The lowest verified loss under the **new target contract** is **48.170707378** at chain **2**. Search budgets ended without an optimizer-convergence certificate. This is an experimental calibration, not an adopted reference.

The empirical wealth/earnings target is **4.45838713455674** (pooled 2005/07 PSID numerator excluding business/farm equity, other real estate, and vehicles; catch-all other assets retained). The numerical objective weight remains **7.595098472533724** as a controlled sensitivity; a new standard error is unavailable. The other 13 target rows and model observers are unchanged. The bequest/wealth denominator compatibility is unresolved.

Complete target fingerprint: `c7a3d185668122e508a6c322bc5ef0715ebb0ecb23948c8d9b184ee25d1cde70`. Complete target-plus-weight fingerprint: `f762ebb5684ab30487b3b8b64fc10977fda396b520035d91c0c5c803255f88e4`. All ten chains share the same authenticated selected source, timing manifest, normalization source pins, and start table. The [verification receipt](verification.json) records each file hash and gate.

## All ten chains

| Chain | Verified new loss | Evaluations | Full GE cases | Search stop |
|---:|---:|---:|---:|---|
| 0 | 52.482365 | 110 | 109 | `native_evaluation_budget_exhausted` |
| 1 | 53.616273 | 108 | 107 | `native_evaluation_budget_exhausted` |
| 2 | 48.170707 | 113 | 112 | `native_evaluation_budget_exhausted` |
| 3 | 54.696044 | 113 | 112 | `native_evaluation_budget_exhausted` |
| 4 | 53.292329 | 109 | 108 | `native_evaluation_budget_exhausted` |
| 5 | 56.376629 | 108 | 107 | `native_evaluation_budget_exhausted` |
| 6 | 55.898691 | 106 | 105 | `native_evaluation_budget_exhausted` |
| 7 | 55.753987 | 106 | 105 | `native_evaluation_budget_exhausted` |
| 8 | 52.583841 | 111 | 110 | `native_evaluation_budget_exhausted` |
| 9 | 54.092429 | 108 | 107 | `native_evaluation_budget_exhausted` |

Winner: [authoritative experimental target table](winner_target_fit.csv), [all 31 parameters](winner_parameters.csv), [selected native report](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fixed_reference_economics_20260928/alternative_wealth_local_20261003_v1/overnight_two_wave/chain_02/native_postcheck/selected_postcheck/phase_b_ge/selected_root), [all-chain machine-readable CSV](all_chains.csv). The native `target_fit.csv` inside the selected report records the **old 6.92658379107299 target** and is diagnostic. `target_fit_new_contract.csv` beside it and the winner table here are authoritative for this experiment. The standard 17 plots diagnose model policies, distributions, and markets; inherited target markers, if any, do not override the new target table.

## Winner target fit: all 14 rows

| Moment | Role | Target | Model | Gap | Weight | Loss contribution |
|---|---|---:|---:|---:|---:|---:|
| initial_normalization | normalization | 2.1 | 2.1 | 1.42031e-08 | — | — |
| cps_childlessness | scored | 0.198279 | 0.20127 | 0.00299115 | 35532.3 | 0.317906 |
| cps_exactly_one | scored | 0.213655 | 0.213233 | -0.00042261 | 26952.8 | 0.00481374 |
| nchs_mean_age | scored | 25.9763 | 25.9532 | -0.0231013 | 139.828 | 0.074622 |
| nchs_share30 | validation | 0.249278 | 0.231581 | -0.0176971 | 0 | 0 |
| wealth_earnings | scored | 4.45839 | 6.56173 | 2.10335 | 7.5951 | 33.6012 |
| bequest_wealth | scored | 0.00729102 | 0.00687019 | -0.00042083 | 5.16529e+06 | 0.914761 |
| old_dispersion | validation | 3.51594 | 2.90173 | -0.61421 | 0 | 0 |
| mean_rooms | scored | 5.72943 | 5.87652 | 0.147084 | 128.021 | 2.76958 |
| ownership_30_55 | scored | 0.67626 | 0.662681 | -0.0135793 | 2339.36 | 0.43137 |
| first_birth_rooms | scored | 1.465 | 1.32951 | -0.135487 | 137.565 | 2.52524 |
| family_rooms | validation | 0.3851 | 0.304636 | -0.080464 | 0 | 0 |
| recent_parent_ownership | scored | 0.127608 | 0.12707 | -0.000538663 | 27055.8 | 0.00785047 |
| early_fertility | scored | 0.809528 | 0.53524 | -0.274288 | 100 | 7.52337 |

## Winner parameters: all 31 rows

| Parameter | Estimate | Lower | Upper | Near bound | Status |
|---|---:|---:|---:|---|---|
| H0 | 6.41684 | 0.2 | 80 | False | derived calibrated housing supply coefficient at N0=1 |
| beta_annual | 0.964495 | 0.94 | 0.99 | False | free in experimental utility calibration |
| chi | 1.04839 | 0.1 | 5 | False | free in experimental utility calibration |
| first_birth_fixed_cost | 0.355747 | 0 | 8 | False | free in experimental utility calibration |
| kappa_fert | 0.122913 | 0.02 | 50 | True | free in experimental utility calibration |
| kappa_fert_continuation | 0.394252 | 0.02 | 50 | True | free in experimental utility calibration |
| theta0 | 0.106175 | 0 | 8 | False | free in experimental utility calibration |
| delta_alpha_jump | 0 | — | — | — | fixed zero under experimental utility contract |
| child_benefit_curvature | 0.106637 | 0 | 0.8 | False | free in experimental utility calibration |
| tenure_choice_kappa | 0.0140426 | 0.001 | 0.1 | False | free in experimental utility calibration |
| psi_child | 0.184318 | 0.01 | 0.5 | False | free in experimental utility calibration |
| child_benefit_CRRA_coefficient | 0.164663 | — | — | — | derived from fixed benefit and proposed curvature |
| theta1 | 0.00819308 | — | — | — | fixed external restriction |
| sigma | 2 | — | — | — | fixed |
| alpha_cons | 0.733 | — | — | — | fixed CEX childless expenditure share |
| delta_alpha | 0 | — | — | — | fixed zero under experimental utility contract |
| h_P | 2.50802 | 0.1 | 2.6 | False | free in experimental utility calibration |
| utility_reference_rent | 0.110466 | — | — | — | retained inactive normalization; compensation off |
| q_annual | 0.02 | — | — | — | author-retained 2% annual real rate |
| financed_share | 0.8 | — | — | — | inherited credit contract |
| housing_supply_elasticity | 0.63 | — | — | — | fixed provisional external mapping |
| payroll_tax | 0.0802807 | — | — | — | derived from adopted pension ratio |
| pension_period | 0.917784 | — | — | — | balanced PAYGO |
| annual_depreciation | 0.0141614 | — | — | — | adopted |
| period_depreciation | 0.0554538 | — | — | — | compounded |
| annual_property_tax | 0.0105984 | — | — | — | adopted |
| period_property_tax | 0.0423934 | — | — | — | linear period convention |
| selling_cost | 0.06 | — | — | — | retained |
| rental_cap | 6 | — | — | — | retained provisional |
| wealth_grid_nodes | 120 | — | — | — | retained exact grid |
| income_states | 9 | — | — | — | retained B15 |

Validation for every chain: 14 experimental target rows; 31 reported parameters, with all ten free estimates inside their bounds; 17 actual PNGs whose hashes match the exact repeat; identical experimental selected/repeat tables; fresh postcheck whose saved search hash matches its search receipt; original/new target arithmetic; search/native loss agreement within 1e-8. See `verification.json` for per-chain paths, hashes, and losses.
