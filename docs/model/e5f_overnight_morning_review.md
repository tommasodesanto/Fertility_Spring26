---
title: "What the overnight checks change"
subtitle: "Housing, tenure, and fertility: discussion brief"
author: "Prepared for Tommaso De Santo"
date: "5 September 2026"
---

# The current assessment

**The supply-versus-credit result survives the saving-optimizer check. Ownership measurement, late-life tenure, and weak parameter identification are now the main issues for discussion.** Production remains September 4 `task_010`, with twelve-moment loss **30.48297**, eleven estimated parameters, housing-supply elasticity **0.63**, and tenure dispersion **0.005**. Nothing has been promoted or re-estimated overnight.

The independent reconstruction reproduces all **75 historical entries exactly**: fifteen fields at each of the five dates from 2007 to 2023. The standard seventeen-figure diagnostic packet now exists for the actual 2023 distribution. The complete production fit and parameter records appear in the appendices.

| Check | New evidence | Implication |
|---|---|---|
| Budget-reporting floors | Excess expenditure affects only \(2.91\times10^{-9}\) of household mass. | Not an economically material explanation of the current calibration miss. |
| Saving optimization | Rare non-global choices exist, but complete global saving changes the tested credit effect by only 0.0000046 percentage points. | This optimizer defect does not explain the weak initial credit-fertility response. |
| Identification | Both twelve- and thirteen-row Jacobians have algebraic rank eleven, but effective rank ten at the declared \(10^{-3}\) relative threshold. | The extra ownership row does not resolve the weak bequest-parameter direction. |
| Young ownership | The reported comparison is 31.10% versus 34.12%; matching annual-age weights and cell boundaries gives 26.98% in one explicit diagnostic. | The three-percentage-point gap is sensitive to measurement, not a cleanly matched age-window fit. |
| Older ownership | Model ownership at ages 80-84 is approximately 99.2%, versus 75.0% in the same ACS household-head sample. | Late-life tenure and housing retention deserve direct validation. |

The first-birth housing response is still the central fit problem: **0.436 rooms versus 0.720 targeted rooms**. The completed sensitivity panel contains a small additional improvement in the original objective, but it barely changes this response. Its best case is reported separately in Appendix B and is not a new production selection.

All new cluster work is complete: the dated-state reconstruction, the 50,000-draw oracle, seven smaller-step historical evaluations, and nine saving-method comparisons. The results below distinguish established findings from checks still needed.

\clearpage

# What global saving changes

The global method checks every interpolation interval when choosing saving, and then solves the complete household problem backward. The experiment holds calibrated parameters and the inherited 2023 population fixed. Each numerical method clears the current housing market separately, with a tighter requested residual of \(10^{-5}\).

| Impact births per household | Production optimizer | Global saving | Difference, pp |
|---|---:|---:|---:|
| Housing supply +20% | +1.290982% | +1.290873% | -0.000109 |
| Dependent-child LTV 95% | +0.017845% | +0.017850% | +0.0000046 |

![Supplemental numerical comparison. Each bar measures the policy change against its own numerical method's baseline; the two panels have different vertical scales.](output/model/e5f_overnight_independent_verification_20260905a/global_saving_receipts/supplemental_global_saving.pdf){width=92%}

The credit effect changes by **0.026% of its own small value**. Tightening market clearing moves the original reported 0.01836% to approximately 0.01785%, a larger change than replacing the saving optimizer. Thus the tiny coefficient should not be given excessive precision, although the supply-versus-credit contrast is stable in this test. Separate clearing reaches identical recorded prices under both methods; the additional fixed-price comparisons therefore coincide at this resolution.

**The remaining numerical priority is discretization.** At baseline, 13.22% of renters are at the six-room cap, including 20.40% of renters with dependent children; 33.49% of owners choose the ten-room maximum. Under supply expansion, the renter shares rise to approximately 26.7% and 40.1%. These are substantial exposures to the housing menu. Adding housing options with option-specific tastes changes the economic choice set, so it must be labeled separately from refining the wealth grid.

All nine cases retain negligible reporting-floor budget excess and valid probabilities. The local credit case has two occupied wealth-monotonicity violations; global saving removes both. Their lower-node mass is only \(1.08\times10^{-7}\) of households. The original-price smoke finds no value decreases across 2.69 million states above the production value cutoff. These are **conditional impact comparisons**, not a revised 2063 transition, a welfare calculation, or a confidence interval.

\clearpage

# Identification after smaller steps

The failed panel postprocessor is now diagnosed. All 23 cases satisfy the original numerical and scientific contracts. Exact queue-flow identities close within \(2.22\times10^{-16}\). A separate comparison between a birth-flow ratio and reconstructed completed fertility differs by at most \(2.62\times10^{-9}\), inside the calibration driver's \(2\times10^{-8}\) measurement tolerance. The ridge planner imposes \(2\times10^{-10}\) on this comparison and rejects eighteen cases. Its code and the production gates have not been weakened; an independent read-only analysis reconstructs the Jacobians from the saved moments.

| Same parameter-coordinate metric | Original 12 rows | Augmented 13 rows |
|---|---:|---:|
| Smallest singular value | 0.560326 | 0.560335 |
| Largest singular value | 579.31 | 853.78 |
| Relative rank at \(10^{-3}\) | 10 of 11 | 10 of 11 |
| Weakest-direction squared loading on \(\theta_1\) | 99.08% | 99.08% |

Here \(\theta_1\) is the wealth shifter in bequest utility. The extra ownership row mainly informs the owner-service premium \(\chi\), already a strongly responsive direction. Adding a row does not destroy information; the larger condition number reflects stronger information in an already strong direction. This is evidence of **weak local identification**, not proof of structural underidentification.

The original panel moves bounded, transformed coordinates by 0.02. That means a positive \(\chi\) step of **8.14%** and a positive supply-intercept step of **12.73%**, not 2% parameter changes. The completed seven-case follow-up uses steps of 0.005 for patience, the owner premium, and the bequest shifter, plus an exactly reproduced anchor.

| Agreement of forward and backward derivatives | Step 0.02 | Step 0.005 |
|---|---:|---:|
| Annual patience \(\beta\) | 0.085 | 0.959 |
| Owner premium \(\chi\) | 0.471 | 0.949 |
| Bequest shifter \(\theta_1\) | 0.972 | 0.877 |

Agreement is the cosine between the two weighted derivative vectors; one means the same direction. Smaller steps improve the first two comparisons substantially, but their central derivative vectors still differ from the larger-step estimates by **35% and 53%**, respectively. This is three checked columns, not a new complete Jacobian or proof of stable rank.

A roughly 2% increase in \(\chi\) raises reported young ownership from **31.10% to 38.97%**, while lowering the first-birth rooms response from **0.4364 to 0.4041**. The mechanism trade-off survives the smaller step. Full target fits, weights, loss contributions, and parameter bounds for all seven cases are retained in `smallstep_receipts/all_target_fits.csv` and `smallstep_receipts/all_parameters.csv`. No calibration proposal was launched.

\clearpage

# Ownership measurement and the demographic handoff

The ACS comparison uses metropolitan household heads, household weights, standard housing structures, and positive room counts, pooled over **2012-2023**. It is not a national 2023 rate. The production model selects complete age nodes 26, 30, and 34 for its 25-34 statistic. The following calculations hold the 2023 model solution fixed and apply the same annual-age ACS weights on both sides.

| Diagnostic measurement | Model | ACS target | Model minus target |
|---|---:|---:|---:|
| Production young statistic, original age aggregation | 31.10% | 34.12% | -3.02 pp |
| Ages 25-34, rates constant over \([a_j,a_j+4)\) | 26.98% | 34.12% | -7.14 pp |
| Ages 25-34, linear interpolation between nodes | 30.22% | 34.12% | -3.90 pp |
| Ages 80-84, rates constant over model cells | 99.23% | 74.97% | +24.26 pp |

These are measurement sensitivities, not newly approved targets: the model does not identify annual within-cell choices. Both aligned young comparisons remain below the data; the old-age discrepancy is large under either interpolation convention.

![Supplemental lifecycle validation. Lines connect model age nodes; the ACS series is pooled across years. The rooms reference is an aggregate target, not an age-specific target.](output/model/e5f_overnight_independent_verification_20260905a/supplemental_figures/lifecycle_validation.pdf){width=95%}

The entry jump is a different issue. Historical household totals and age shares are imposed, while future entrants arrive from an inherited birth queue. In the baseline,

\[
E_{2023}=0.0143371,\qquad E_{2027}=0.0639084=0+1\times B_{2023}.
\]

Thus first-period entry rises by a factor of **4.46** when the age-raking adjustment stops. The increase is entirely the handoff to the due historical queue, not a response to policy births in 2023. A 2023 birth affects entry only twenty years later in this contract.

Report these as **conditional household mechanisms from a fitted 2023 state**. A resident-population forecast requires a justified household-formation and entry transition.

\clearpage

# What to decide when we meet

1. **Set one age-measurement convention before activating young ownership.** Keep the current twelve rows during that decision. Adding a well-defined holdout is useful, but the household-head correction alone does not finish alignment of age, geography, and time.

2. **Use older ownership and housing retention to assess the bequest block.** Nearly universal late-life ownership is a clear validation issue. Fixing \(\theta_1\) at its current value removes the weak numerical direction in a column-deletion calculation, but that calculation supplies no external economic justification for the restriction. We should choose evidence or a disclosed restriction, rather than silently dropping a parameter or target.

3. **Keep the first-birth room response at the center of calibration diagnosis.** The new panel's best twelve-row score improves from 30.483 to 29.386, but the room response rises only from 0.4364 to 0.4394. A smaller objective is not resolution of the mechanism's main empirical miss.

4. **Separate initial policy mechanisms from the household-entry transition.** The existing supply and credit contrast is conditional on this model and starting distribution. Long-horizon aggregate births additionally reflect the queue and household counts. The rebate result also depends on the fiscal comparison: rebated 2% versus rebated 1%, not the unrebated status quo.

5. **Move the next numerical check to wealth grids, transaction interpolation, and the housing menu.** Global saving has negligible effects on the two initial policy responses tested. The substantial mass at rental and owner size limits makes housing-menu sensitivity more relevant than another search at the current discretization. Keep numerical refinement distinct from adding economic housing choices.

The author-owned manuscript remains untouched. No new target, weight, parameter restriction, or production benchmark has been adopted. The H128 baseline remains unaccepted; no unchanged long-horizon continuation was launched.

**Evidence and reproduction.** The complete working packet is in `output/model/e5f_overnight_independent_verification_20260905a/`. The numerical checks run against the frozen production code, with all fifty recorded source hashes verified. The full audit remains the detailed reference; its overnight addendum records changes to previously pending findings.

The [panel receipts](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_overnight_independent_verification_20260905a/panel_receipts/summary.json), [numerical receipts](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_overnight_independent_verification_20260905a/numerical_full/summary.json), [annual-age measurement receipts](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_overnight_independent_verification_20260905a/panel_receipts/ownership_annual_age_receipts.csv), and [run plan](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/docs/model/e5f_overnight_verification_plan_20260905.md) retain exact inputs, assumptions, and stop rules. Complete fit and parameter tables follow.


\clearpage

# Appendix A. Complete production target fit

The twelve-row objective is **30.4829667**. Gaps equal model minus target; each contribution is weight times squared gap. The young-ownership row is additional and default-off.

| Moment | Target | Model | Gap | Weight | Loss |
|---|---:|---:|---:|---:|---:|
| Completed fertility | 1.918000 | 1.922907 | 0.004907 | 1425.74 | 0.0343 |
| Childlessness | 0.188000 | 0.189407 | 0.001407 | 17180.7 | 0.0340 |
| Age at first birth | 26.044627 | 26.256032 | 0.211404 | 44.4444 | 1.9863 |
| First births age 30+ | 0.260327 | 0.237579 | -0.022748 | 10000 | 5.1748 |
| First-birth rooms response | 0.720246 | 0.436418 | -0.283828 | 137.565 | 11.0821 |
| Rooms, 3+ minus 1-2 children | 0.367700 | 0.403441 | 0.035741 | 2958.51 | 3.7793 |
| Family ownership gap | 0.167662 | 0.161699 | -0.005963 | 14229.6 | 0.5060 |
| Ownership rate, ages 30-55 | 0.575472 | 0.544488 | -0.030984 | 1207.85 | 1.1596 |
| Mean occupied rooms | 5.779970 | 6.317291 | 0.537321 | 11.9732 | 3.4568 |
| Wealth / annual earnings | 6.873100 | 6.932652 | 0.059552 | 6.28767 | 0.0223 |
| Bequest flow / wealth | 0.008800 | 0.008433 | -0.000367 | 5.16529e+06 | 0.6971 |
| Older wealth p90 / p50 | 3.448111 | 3.236511 | -0.211600 | 56.9598 | 2.5503 |

At this exact parameter vector, young ownership is 0.31098385 against 0.34116609. Its synthetic weight is 3436.594375 and its additional loss is 3.1306261, giving a thirteen-row score of 33.6135928. This is not a deterioration in the twelve-row fit.

# Complete parameter comparison

All eleven free parameters and both external restrictions are shown. The production and panel-003 cases differ only in annual patience; the old fertility-preference intercept is re-normalized to replacement in each case.

| Parameter | Production | Panel 003 | Bounds / restriction | Near bound? |
|---|---:|---:|---|---|
| Annual $\beta$ | 0.995117 | 0.995739 | [0.94, 0.9995] | No |
| First-birth $\kappa_1$ | 2.168173 | 2.168173 | [0.02, 50] | No |
| Later-birth $\kappa_C$ | 1.736471 | 1.736471 | [0.02, 50] | No |
| Owner premium $\chi$ | 1.043472 | 1.043472 | [0.1, 5] | No |
| Supply intercept $H_0$ | 14.562959 | 14.562959 | [0.2, 80] | No |
| Bequest strength $\theta_0$ | 0.528428 | 0.528428 | [0, 8] | No |
| Bequest shifter $\theta_1$ | 0.107249 | 0.107249 | [0.02, 16] | Raw scale only |
| Per-child room floor | 0.282210 | 0.282210 | [0.1, 1.8] | No |
| First-birth fixed cost | 4.559138 | 4.559138 | [0, 8] | No |
| First-child room jump | 0.364931 | 0.364931 | [0, 0.5] | No |
| Child-value change | -0.328714 | -0.328714 | [-1.5, 0.2] | No |
| Tenure $\kappa$ | 0.005000 | 0.005000 | Externally fixed | Not estimated |
| Supply elasticity $\eta$ | 0.630000 | 0.630000 | Externally fixed | Not estimated |

The stored near-bound flag for $\theta_1$ uses raw units. Its position is 25.1% of the way through the logarithmic search interval, so it is not mechanically stuck at the lower search edge. Production $\psi_{2007}=0.28801685$ and $\psi_{2023}=-0.04069705$; panel 003 uses 0.28667671 and -0.04203720.

\clearpage

# Appendix B. Complete unpromoted panel-003 fit

This is the best observed case in the original 23-case young-ownership panel, not a newly launched calibration. Its original twelve-row score is **29.3858837** and augmented score is **32.3707373**. That is a 3.60% improvement in the original objective relative to production; the first-birth room shortfall is still 0.28090 rooms. All parameters, bounds, and restrictions appear in Appendix A.

| Moment | Target | Model | Gap | Weight | Loss |
|---|---:|---:|---:|---:|---:|
| Completed fertility | 1.918000 | 1.921701 | 0.003701 | 1425.74 | 0.0195 |
| Childlessness | 0.188000 | 0.189212 | 0.001212 | 17180.7 | 0.0252 |
| Age at first birth | 26.044627 | 26.247452 | 0.202825 | 44.4444 | 1.8284 |
| First births age 30+ | 0.260327 | 0.236974 | -0.023354 | 10000 | 5.4539 |
| First-birth rooms response | 0.720246 | 0.439350 | -0.280896 | 137.565 | 10.8542 |
| Rooms, 3+ minus 1-2 children | 0.367700 | 0.404458 | 0.036759 | 2958.51 | 3.9976 |
| Family ownership gap | 0.167662 | 0.164819 | -0.002843 | 14229.6 | 0.1150 |
| Ownership rate, ages 30-55 | 0.575472 | 0.546824 | -0.028648 | 1207.85 | 0.9913 |
| Mean occupied rooms | 5.779970 | 6.320025 | 0.540054 | 11.9732 | 3.4921 |
| Wealth / annual earnings | 6.873100 | 7.015451 | 0.142351 | 6.28767 | 0.1274 |
| Bequest flow / wealth | 0.008800 | 0.008417 | -0.000383 | 5.16529e+06 | 0.7573 |
| Older wealth p90 / p50 | 3.448111 | 3.274138 | -0.173972 | 56.9598 | 1.7240 |
| Young ownership (diagnostic) | 0.341166 | 0.311695 | -0.029471 | 3436.59 | 2.9849 |

The thirteenth-row scale is **5% of its target**, not an empirical standard error. These objective values are diagnostic scores, not specification-test statistics. Source and target fingerprints differ from the production record only as explicitly documented for the additional default-off row and its code profile. The household numerical kernels are unchanged in these inherited panel evaluations.

The working packet includes the complete 23-case tables:

- `panel_receipts/all_target_fits.csv`
- `panel_receipts/all_parameters.csv`

No case has been promoted, and no ridge proposal has been submitted.
