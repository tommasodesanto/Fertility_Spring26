# Numerical pair: completed diagnostic

Torch 18753562 completed both arms successfully. The same proposed parameters were evaluated in the original one-birth model. All native scientific gates and pre-specified paired screens passed. The reference remains **2007 stationary reference — block0506, September 28 verified export**; no candidate is promoted.

The joint parameter step reduced loss from 19.581 to 13.774 / 13.779, a 29.65% reduction. The derivative prediction was 13.335. Children by 25 remain about 0.533 against 0.810. Most of the remaining scored loss is now the early-fertility row.

## Numerical comparison

| Starting guess | Stationary solves | Solving time, minutes | Total time, minutes | Loss | Final psi | Completed fertility |
|---|---:|---:|---:|---:|---:|---:|
| Retained | 7 | 18.431 | 19.541 | 13.774 | 0.131 | 2.100 |
| Predicted | 3 | 8.558 | 9.660 | 13.779 | 0.131 | 2.100 |

The predicted start saved 50.56% of total time and 53.57% of solving time at this point. It uses the unchanged fertility tolerance of 5.000e-04. The retained arm happened to finish closer to exactly 2.1. Loss differs by 4.388e-03; maximum square-root-weight-scaled scored-moment difference is 7.216e-03, below the preset 0.01 limit. Price, psi and every validation screen also pass. This is a one-point numerical comparison, not a statistical equivalence test or a general speed benchmark.

## Complete fit

The target, retained-start model, gap, weight and loss columns report the full active contract. The predicted-start model is shown beside it; its full gap and loss columns remain in the linked raw table.

Early fertility counts children ever born capped at three. The family-rooms
validation row is a cross-sectional comparison by resident-child count, using
the retained model dependent-child proxy; it is not a second-birth event study.
All inherited cohort/period, age-cell and observer approximations remain in place.

### Targeted rows

| Moment | Target | Reference | Retained-start model | Predicted-start model | Retained gap | Weight | Retained loss |
|---|---:|---:|---:|---:|---:|---:|---:|
| Childlessness, ages 40–44 | 0.198 | 0.201 | 0.200 | 0.200 | 1.966e-03 | 35532.3 | 0.137 |
| Exactly one child among mothers, ages 40–44 | 0.214 | 0.209 | 0.211 | 0.211 | -3.128e-03 | 26952.82 | 0.264 |
| Mean first-birth age | 25.976 | 25.933 | 25.951 | 25.951 | -0.025 | 139.828 | 0.089 |
| Wealth / earnings | 6.927 | 6.326 | 6.585 | 6.585 | -0.342 | 7.595 | 0.888 |
| Bequests / wealth | 7.291e-03 | 7.057e-03 | 7.089e-03 | 7.089e-03 | -2.024e-04 | 5165289 | 0.212 |
| Mean rooms | 5.729 | 5.848 | 5.816 | 5.816 | 0.087 | 128.021 | 0.961 |
| Ownership, ages 30–55 | 0.676 | 0.655 | 0.661 | 0.661 | -0.015 | 2339.362 | 0.510 |
| First-birth rooms response | 1.465 | 1.622 | 1.576 | 1.576 | 0.111 | 137.565 | 1.694 |
| Recent-parent ownership gap | 0.128 | 0.120 | 0.120 | 0.120 | -7.116e-03 | 27055.82 | 1.370 |
| Children by age 25 | 0.810 | 0.535 | 0.533 | 0.533 | -0.277 | 100 | 7.649 |

### Untargeted validation rows

| Moment | Target | Reference | Retained-start model | Predicted-start model | Retained gap | Weight | Retained loss |
|---|---:|---:|---:|---:|---:|---:|---:|
| First births at 30+ | 0.249 | 0.224 | 0.224 | 0.224 | -0.025 | 0 | 0 |
| Wealth/income p90/p50, ages 76–84 | 3.516 | 3.069 | 2.999 | 2.999 | -0.517 | 0 | 0 |
| Rooms: 3+ versus 1–2 children, ages 30–55 | 0.385 | 0.353 | 0.391 | 0.391 | 5.784e-03 | 0 | 0 |

### Separate replacement normalization

| Moment | Target | Reference | Retained-start model | Predicted-start model | Retained gap | Weight | Retained loss |
|---|---:|---:|---:|---:|---:|---:|---:|
| Completed fertility | 2.1 | 2.100 | 2.100 | 2.100 | -1.662889e-06 | — | — |

## Parameters and bounds

The ten proposed calibration coordinates are identical across arms. No bounds changed. Near-bound means within 1% of the full inherited interval, which is wide for the two fertility scales; it does not mean either parameter is numerically equal to its bound.

| Parameter | Estimate | Lower | Upper | Near bound |
|---|---:|---:|---:|---|
| H0 | 6.224 | 0.2 | 80 | False |
| beta_annual | 0.966 | 0.94 | 0.99 | False |
| chi | 1.094 | 0.1 | 5 | False |
| first_birth_fixed_cost | 0.534 | 0 | 8 | False |
| kappa_fert | 0.155 | 0.02 | 50 | True |
| kappa_fert_continuation | 0.293 | 0.02 | 50 | True |
| theta0 | 0.127 | 0 | 8 | False |
| delta_alpha_jump | 0.127 | 0 | 0.25 | False |
| child_benefit_curvature | 0.062 | 0 | 0.8 | False |
| tenure_choice_kappa | 0.012 | 0.001 | 0.1 | False |

Child benefit is separately normalized and differs slightly across arms, as shown above. All remaining external restrictions and derived values appear in both complete 31-row parameter tables.

- [Retained-start full fit](run_v1/retained_start/case/target_fit.csv) and [all parameters/restrictions](run_v1/retained_start/case/parameters.csv).
- [Predicted-start full fit](run_v1/predicted_start/case/target_fit.csv) and [all parameters/restrictions](run_v1/predicted_start/case/parameters.csv).
- [All predictions versus actual outcomes](run_v1/predictions_vs_actual.csv) and [paired numerical screens](run_v1/paired_comparison.json).
- [Dispatch receipt, including all ten completed stationary solves](run_v1/dispatch_receipt.json).

## Interpretation and limits

The observed local joint step improves wealth, ownership and housing-response fits, while early fertility falls slightly. The older-age dispersion validation moment worsens from 3.069 to 2.999 against 3.516; this zero-weight row is retained. The two additional validation rows improve. All scientific gates pass, but this does not establish precise identification, global optimality, or suitability for the 2023 transition.

Both cases retain all 17 standard plots. The predicted-start set was visually inspected: the existing high-wealth age-30 housing drop, small ownership downturn and retirement wealth kink remain; the market residual is 9.86e-7. The retained-start plots are collected, with native/controller checks, but not claimed as a separately completed visual review. Buyer-conditional policy-panel scope remains unchanged. Checkpoints remain on Torch and outside Git.

This successful pair makes the joint step and initial-psi prediction useful candidates for the next bounded calibration experiment. No final exact repeats or further search are part of this budget. Derivative-based recalibration of the changed two-birth household problem would require fresh derivatives; its currently running fixed-parameter test does not. The other chat continues to use its frozen reference.

The inherited filename `lifecycle_2023.csv` in each native diagnostic packet contains this stationary solution's profile. It is not a 2023 transition result; no transition was run. Standard filenames were retained.

The author's supplemental four-panel lifecycle comparison is regenerated by
`build_lifecycle_fit.py` on Torch from the retained-start checkpoint, with no
model solves. It preserves the earlier children, ownership, rooms and net-worth
definitions and empirical samples from `overnight_calibration_20260928/morning_fits/`.
`lifecycle_fit.png` matches the single model column in the latest chat table;
`model_profiles.csv` and `lifecycle_fit_qa.json` retain its numbers and source
checks. These descriptive age profiles are untargeted, apart from the separately
marked exact-age-25 fertility moment. The 17 standard diagnostics remain intact.

Torch 18756525 completed this saved-state render in 34 seconds. Checkpoint and
receipt pins, empirical-table identity, exact age-25 fertility and aggregate
wealth/earnings replay passed; the collected four-panel PNG was visually checked.
