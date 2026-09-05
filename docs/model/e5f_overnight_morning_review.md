---
title: "What the overnight checks change"
subtitle: "Housing, tenure, and fertility: discussion brief"
author: "Prepared for Tommaso De Santo"
date: "5 September 2026"
---

# The current assessment

**The supply-versus-credit result survives the saving-optimizer and wealth-grid checks. Ownership measurement, housing opportunities, and weak parameter identification are now the main issues for discussion.** Production remains September 4 `task_010`, with twelve-moment loss **30.48297**, eleven estimated parameters, housing-supply elasticity **0.63**, and tenure dispersion **0.005**. Nothing has been promoted or re-estimated overnight.

The independent reconstruction reproduces all **75 historical entries exactly**: fifteen fields at each of the five dates from 2007 to 2023. The standard seventeen-figure diagnostic packet now exists for the actual 2023 distribution. The complete production fit and parameter records appear in the appendices.

| Check | New evidence | Implication |
|---|---|---|
| Budget-reporting floors | Excess expenditure affects only \(2.91\times10^{-9}\) of household mass. | Not an economically material explanation of the current calibration miss. |
| Saving optimization | Rare non-global choices exist, but complete global saving changes the tested credit effect by only 0.0000046 percentage points. | This optimizer defect does not explain the weak initial credit-fertility response. |
| Wealth-grid resolution | Moving from 120 to 239 points changes the supply birth effect by -0.19% and the credit effect by +3.61% of their respective values. | The contrast survives; tiny credit effects need matched baselines and restrained precision. |
| Identification | Both twelve- and thirteen-row Jacobians have algebraic rank eleven, but effective rank ten at the declared \(10^{-3}\) relative threshold. | The extra ownership row does not resolve the weak bequest-parameter direction. |
| Young ownership | The reported comparison is 31.10% versus 34.12%; matching annual-age weights and cell boundaries gives 26.98% in one explicit diagnostic. | The three-percentage-point gap is sensitive to measurement, not a cleanly matched age-window fit. |
| Older ownership | Model ownership at ages 80-84 is approximately 99.2%, versus 75.0% in the same ACS household-head sample. | Late-life tenure and housing retention deserve direct validation. |

The first-birth housing response is still the central fit problem: **0.436 rooms versus 0.720 targeted rooms**. The completed sensitivity panel contains a small additional improvement in the original objective, but it barely changes this response. Its best case is reported separately in Appendix B and is not a new production selection.

All planned cluster work is complete, including the nested-grid check and separate rental-cap sensitivity. Larger rentals leave supply's birth response near **1.29%** but reverse its ownership response. The remaining decisions concern measurement, identification, and economic restrictions.

\clearpage

# What global saving changes

The global method checks every interpolation interval when choosing saving, and then solves the complete household problem backward. The experiment holds calibrated parameters and the inherited 2023 population fixed. Each numerical method clears the current housing market separately, with a tighter requested residual of \(10^{-5}\).

| Impact births per household | Production optimizer | Global saving | Difference, pp |
|---|---:|---:|---:|
| Housing supply +20% | +1.290982% | +1.290873% | -0.000109 |
| Dependent-child LTV 95% | +0.017845% | +0.017850% | +0.0000046 |

![Supplemental numerical comparison. Each bar measures the policy change against its own numerical method's baseline; the two panels have different vertical scales.](output/model/e5f_overnight_independent_verification_20260905a/global_saving_receipts/supplemental_global_saving.pdf){width=92%}

The credit effect changes by **0.026% of its own small value**. Tightening market clearing moves the original reported 0.01836% to approximately 0.01785%, a larger change than replacing the saving optimizer. Thus the tiny coefficient should not be given excessive precision, although the supply-versus-credit contrast is stable in this test. Separate clearing reaches identical recorded prices under both methods; the additional fixed-price comparisons therefore coincide at this resolution.

**The housing-size limits are economically consequential restrictions.** At baseline, 13.22% of renters are at the six-room cap, including 20.40% of renters with dependent children; 33.49% of owners choose the ten-room maximum. Under supply expansion, the renter shares rise to approximately 26.7% and 40.1%. These are substantial exposures to the housing menu. Adding housing options with option-specific tastes changes the economic choice set, so it must be labeled separately from refining the wealth grid.

All nine cases retain negligible reporting-floor budget excess and valid probabilities. The local credit case has two occupied wealth-monotonicity violations; global saving removes both. Their lower-node mass is only \(1.08\times10^{-7}\) of households. The striking first-birth spike at negative wealth in the age-42 renter plot has **exactly zero household mass** in all three cases; all local/global first-birth probability differences above 0.1 are also at unoccupied states. The original-price smoke finds no value decreases across 2.69 million states above the production value cutoff. These are **conditional impact comparisons**, not a revised 2063 transition, a welfare calculation, or a confidence interval.

\clearpage

# Wealth-grid refinement

The finer grid retains all 120 original wealth points and adds their interval
midpoints, giving **239 points**. Initial household mass stays exactly at the
original points. Both resolutions use global saving and separately clear each
current policy market. Parameters, inherited history, wealth endpoints, and
housing choices are unchanged.

| Impact births per household | 120 points | 239 points | Change in effect |
|---|---:|---:|---:|
| Housing supply +20% | +1.290873% | +1.288387% | -0.193% |
| Dependent-child LTV 95% | +0.017850% | +0.018494% | +3.61% |

The last column is the relative change in each policy effect, not its difference
in percentage points. **Supply still moves births much more than credit.** The
credit coefficient remains small at both resolutions, though its last digits
are less stable. All three fine-grid markets pass, with maximum residual
\(2.29\times10^{-6}\), valid probabilities, preserved population, and no occupied
wealth-value declines above the declared screen.

![Supplemental grid comparison. Panels use different vertical scales; every policy effect is measured against its own grid's baseline.](output/model/e5f_overnight_independent_verification_20260905a/continuation_receipts/supplemental_nested_grid.pdf){width=92%}

Refinement changes baseline births per household by **0.01748%** and baseline
ownership by **0.00443 pp**. The baseline birth-level change is similar in size
to the credit effect. Comparing a fine-grid policy with the coarse-grid baseline
would therefore give a misleading answer; the table uses matched baselines.

None of the predeclared materiality thresholds is crossed, so this branch stops
at 239 points. This is a bounded resolution check, not a convergence certificate.
The upper wealth bound of 30, the historical path, the matched 2019-2023
first-birth rooms target, and long-run policy responses have not been retested.
Changing rental or ownership housing opportunities is a separate economic
sensitivity, described separately from this numerical check.

\clearpage

# Rental opportunities: an economic sensitivity

The model caps rentals at six rooms. The existing ACS/MMS housing-stock packet
places **9.17% of renter households in seven-or-more-room units**. Saved weighted
cells reproduce that share in a descriptive sample that differs from the active
ownership target. It motivates validation without identifying a replacement limit.

The eight-room experiment fixes all eleven baseline estimates, owner choices,
external restrictions, the 239-point wealth grid, and the inherited 2023 household
distribution. Each policy uses its own rental-limit baseline. Eight is an
**uncalibrated economic choice**, not a numerical repair or production default.

| Impact births per household | Six-room limit | Eight-room limit |
|---|---:|---:|
| Supply +20% | +1.288387% | +1.292605% |
| Dependent-child LTV 95% | +0.018494% | +0.013773% |

![Supplemental rental-cap sensitivity. Each policy uses its own baseline; panels have different units and vertical scales.](output/model/e5f_overnight_independent_verification_20260905a/continuation_receipts/supplemental_rental_cap.pdf){width=92%}

**Supply's birth response barely changes, while its ownership effect moves from
+1.770 pp to -0.232 pp.** Credit's birth effect is **25.5% smaller** and remains
very small; ownership still rises by **0.626 pp**. The contrast between the birth
effects survives, but the exact credit coefficient and tenure predictions depend
on the housing opportunity set.


At baseline, allowing eight rooms lowers ownership from **63.57% to 61.60%**,
increases rooms per household by **0.386%**, and changes births by **-0.0198%**.
Prices and tenure composition adjust together. These changes show sensitivity
to the housing opportunity set, without establishing a better joint fit.

The matched 2019-2023 first-birth room response, full calibration score, history,
long-run transitions, and welfare have not been recomputed. The next decision
is how tenure-specific size distributions and transitions should discipline
the model.

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

The ACS comparison uses metropolitan household heads, household weights, standard housing structures, and positive room counts, pooled over **2012-2023**. It is not a national 2023 rate. The historical age bridge assigns annual ages 26-29 to node 26, 30-33 to node 30, and 34-37 to node 34. The model selects these three complete nodes for its reported 25-34 statistic. The following calculations hold the 2023 solution fixed; the aligned rows use the same annual-age ACS weights on both sides.

| Diagnostic measurement | Model | ACS target | Model minus target |
|---|---:|---:|---:|
| Production young statistic, original age aggregation | 31.10% | 34.12% | -3.02 pp |
| Full model cells, ages 26-37, common ACS age weights | 31.18% | 39.51% | -8.33 pp |
| Ages 25-34, rates constant over \([a_j,a_j+4)\) | 26.98% | 34.12% | -7.14 pp |
| Ages 25-34, linear interpolation between nodes | 30.22% | 34.12% | -3.90 pp |
| Ages 80-84, rates constant over model cells | 99.23% | 74.97% | +24.26 pp |

Moving the ACS window from 25-34 to the full model cells, 26-37, increases the empirical rate by **5.39 pp**; changing model age weights then adds only **0.08 pp** to model ownership. The age window drives this difference. These are diagnostic measurements, not newly approved targets: annual within-cell choices are not identified. The old-age discrepancy remains large under either convention.

![Supplemental lifecycle validation. Lines connect model age nodes; the ACS series is pooled across years. The rooms reference is an aggregate target, not an age-specific target.](output/model/e5f_overnight_independent_verification_20260905a/supplemental_figures/lifecycle_validation.pdf){width=95%}

The entry jump is a different issue. Historical household totals and age shares are imposed, while future entrants arrive from an inherited birth queue. In the baseline,

\[
E_{2023}=0.0143371,\qquad E_{2027}=0.0639084=0+1\times B_{2023}.
\]

Thus first-period entry rises by a factor of **4.46** when the age-raking adjustment stops. The increase is entirely the handoff to the due historical queue, not a response to policy births in 2023. A 2023 birth affects entry only twenty years later in this contract.

Report these as **conditional household mechanisms from a fitted 2023 state**. A resident-population forecast requires a justified household-formation and entry transition.

\clearpage

# What to decide when we meet

1. **Set one age-measurement convention before activating young ownership.** Review the existing age-bounded targets too: the statistic labeled 30-55 covers full cells through 57. A literal 30-55 comparison with common ACS weights gives model ownership 53.47% versus 57.55%, instead of the stored model rate 54.45%. Keep the current twelve rows during this decision; age, geography, and time need one explicit contract.

2. **Use older ownership and housing retention to assess the bequest block.** Nearly universal late-life ownership is a clear validation issue. Fixing \(\theta_1\) at its current value removes the weak numerical direction in a column-deletion calculation, but that calculation supplies no external economic justification for the restriction. We should choose evidence or a disclosed restriction, rather than silently dropping a parameter or target.

3. **Keep the first-birth room response at the center of calibration diagnosis.** The new panel's best twelve-row score improves from 30.483 to 29.386, but the room response rises only from 0.4364 to 0.4394. A smaller objective is not resolution of the mechanism's main empirical miss.

4. **Separate initial policy mechanisms from the household-entry transition.** The existing supply and credit contrast is conditional on this model and starting distribution. Long-horizon aggregate births additionally reflect the queue and household counts. The rebate result also depends on the fiscal comparison: rebated 2% versus rebated 1%, not the unrebated status quo.

5. **Use the completed numerical checks to focus validation on housing opportunities.** Global saving and the 120-to-239-point grid check preserve the supply-versus-credit contrast. The rental-cap experiment is an uncalibrated economic sensitivity. Decide which tenure-specific size distributions and transitions should discipline that opportunity set; the wealth ceiling, transaction interpolation, and reconstructed history also remain to be checked.

The author-owned manuscript remains untouched. No new target, weight, parameter restriction, or production benchmark has been adopted. The H128 baseline remains unaccepted; no unchanged long-horizon continuation was launched.

**Evidence and reproduction.** The complete working packet is in `output/model/e5f_overnight_independent_verification_20260905a/`. The numerical checks run against the frozen production code, with all fifty recorded source hashes verified. The full audit remains the detailed reference; its overnight sections record changes to previously pending findings. The `continuation_receipts/` folder contains the age-window decomposition, occupied-state exposure, grid comparison, and rental-cap comparison.

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
| Older wealth/income p90 / p50 | 3.448111 | 3.236511 | -0.211600 | 56.9598 | 2.5503 |

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
| Older wealth/income p90 / p50 | 3.448111 | 3.274138 | -0.173972 | 56.9598 | 1.7240 |
| Young ownership (diagnostic) | 0.341166 | 0.311695 | -0.029471 | 3436.59 | 2.9849 |

The thirteenth-row scale is **5% of its target**, not an empirical standard error. These objective values are diagnostic scores, not specification-test statistics. Source and target fingerprints differ from the production record only as explicitly documented for the additional default-off row and its code profile. The household numerical kernels are unchanged in these inherited panel evaluations.

The working packet includes the complete 23-case tables:

- `panel_receipts/all_target_fits.csv`
- `panel_receipts/all_parameters.csv`

No case has been promoted, and no ridge proposal has been submitted.
