---
title: "What the bounded calibration refinement achieved"
subtitle: "Same model, targets, weights, and external restrictions"
author: "Prepared for Tommaso De Santo"
date: "5 September 2026"
---

# The result for discussion

The best verified candidate has loss **23.153445**, compared with **30.482967**
for the retained production calibration: a **24.04% reduction** in the
same twelve-row objective. Two fresh executions reproduce every fit row,
parameter, preference normalization, and saved numerical history entry exactly.
This is a bounded local improvement, with weak parameter identification still
unresolved. It is recorded for review; production has not been silently replaced.

The first-birth housing response moves from **0.4364**
to **0.4660 rooms**, against **0.7202** targeted.
Mean occupied rooms move from **6.3173**
to **6.3286**, against **5.7800** targeted.
The first-birth response explains **64.7%**
of its target. A lower total objective must be assessed alongside this remaining gap.

The largest gain comes from the rooms contrast between families with three or more
children and families with one or two: its model value changes from
**0.4034** to
**0.3792 rooms**, against
**0.3677** targeted.
Mean occupied rooms remain too high and become slightly higher at this candidate.

| Evaluation stage | Completed cases | Best loss in stage |
|---|---:|---:|
| Inherited-candidate exact repeats | 2 | 29.385884 |
| Small coordinate moves | 23 | 28.587975 |
| Joint parameter moves | 12 | 23.153445 |
| Final exact repeats | 2 | 23.153445 |

All **39** completed evaluations retain the original source and target
fingerprints, twelve moments, eleven free parameters, the 120-point wealth grid,
dated housing-supply elasticity 0.63, and tenure-choice dispersion 0.005. They reconstruct
the old state and full 2007-2023 path. Every case has a terminal checkpoint,
complete target and parameter tables, and the unchanged seventeen-figure packet.
No policy path or long-run response was recalculated in this exercise.

\clearpage

# Where the fit changes

![Supplemental comparison of every objective contribution. Both series use identical target values and weights.](output/model/e5f_morning_refinement_20260905a/report/loss_contributions.pdf){width=100%}

The first round tests small separate parameter moves. The second combines
successful moves and tests a child-space rebalancing: increase the first-child
intercept while reducing the per-child floor. With
$\bar h(m)=\bar h_J+m\bar h_C$ for children at home $m>0$, the proposed changes
$d\bar h_J=d$ and $d\bar h_C=-d/3$ preserve the three-child floor and raise the
first-child floor by $2d/3$. Actual occupied rooms still depend on prices, saving,
tenure choice, and fertility selection; the floor identity is not an identity
for observed housing use.

All ranking uses the complete unchanged objective. No target was removed,
reweighted, or added. The earlier ridge planner's rejection is retained; these
direct candidate evaluations do not invert its unstable Jacobian or establish
precise estimates of the weakly informed bequest parameter.

\clearpage

# Complete target fit: verified candidate

Gaps equal model minus target. Contributions are weight times squared gap.
The weights retain their existing empirical or synthetic interpretations;
this loss is not a specification-test statistic.

| Moment | Target | Model | Gap | Weight | Loss |
|---|---:|---:|---:|---:|---:|
| Completed fertility | 1.918000 | 1.921698 | +0.003698 | 1425.74 | 0.0195 |
| Childlessness | 0.188000 | 0.189859 | +0.001859 | 17180.7 | 0.0594 |
| Mean age at first birth | 26.044627 | 26.251703 | +0.207076 | 44.4444 | 1.9058 |
| First births at age 30+ | 0.260327 | 0.237163 | -0.023165 | 10000 | 5.3660 |
| Rooms response to first birth | 0.720246 | 0.466044 | -0.254202 | 137.565 | 8.8893 |
| Rooms, 3+ minus 1-2 children | 0.367700 | 0.379198 | +0.011499 | 2958.51 | 0.3912 |
| Family ownership gap | 0.167662 | 0.169578 | +0.001917 | 14229.6 | 0.0523 |
| Ownership, ages 30-55 | 0.575472 | 0.546618 | -0.028854 | 1207.85 | 1.0056 |
| Mean occupied rooms | 5.779970 | 6.328594 | +0.548624 | 11.9732 | 3.6038 |
| Wealth / annual earnings | 6.873100 | 7.031196 | +0.158096 | 6.28767 | 0.1572 |
| Bequest flow / wealth | 0.008800 | 0.008485 | -0.000315 | 5.16529e+06 | 0.5139 |
| Older wealth/income p90 / p50 | 3.448111 | 3.303595 | -0.144516 | 56.9598 | 1.1896 |

# Complete target fit: retained production

| Moment | Target | Model | Gap | Weight | Loss |
|---|---:|---:|---:|---:|---:|
| Completed fertility | 1.918000 | 1.922907 | +0.004907 | 1425.74 | 0.0343 |
| Childlessness | 0.188000 | 0.189407 | +0.001407 | 17180.7 | 0.0340 |
| Mean age at first birth | 26.044627 | 26.256032 | +0.211404 | 44.4444 | 1.9863 |
| First births at age 30+ | 0.260327 | 0.237579 | -0.022748 | 10000 | 5.1748 |
| Rooms response to first birth | 0.720246 | 0.436418 | -0.283828 | 137.565 | 11.0821 |
| Rooms, 3+ minus 1-2 children | 0.367700 | 0.403441 | +0.035741 | 2958.51 | 3.7793 |
| Family ownership gap | 0.167662 | 0.161699 | -0.005963 | 14229.6 | 0.5060 |
| Ownership, ages 30-55 | 0.575472 | 0.544488 | -0.030984 | 1207.85 | 1.1596 |
| Mean occupied rooms | 5.779970 | 6.317291 | +0.537321 | 11.9732 | 3.4568 |
| Wealth / annual earnings | 6.873100 | 6.932652 | +0.059552 | 6.28767 | 0.0223 |
| Bequest flow / wealth | 0.008800 | 0.008433 | -0.000367 | 5.16529e+06 | 0.6971 |
| Older wealth/income p90 / p50 | 3.448111 | 3.236511 | -0.211600 | 56.9598 | 2.5503 |

\clearpage

# Complete parameter comparison

| Parameter | Production | New candidate | Bounds / restriction | Near bound? |
|---|---:|---:|---|---|
| Annual patience $\beta$ | 0.995117 | 0.995739 | [0.94, 0.9995] | No |
| First-birth dispersion $\kappa_1$ | 2.168173 | 2.168173 | [0.02, 50] | No |
| Later-birth dispersion $\kappa_C$ | 1.736471 | 1.736471 | [0.02, 50] | No |
| Owner-service premium $\chi$ | 1.043472 | 1.043472 | [0.1, 5] | No |
| Supply intercept $H_0$ | 14.562959 | 14.562959 | [0.2, 80] | No |
| Bequest strength $\theta_0$ | 0.528428 | 0.549189 | [0, 8] | No |
| Bequest shifter $\theta_1$ | 0.107249 | 0.107249 | [0.02, 16] | Raw scale only |
| Per-child room floor | 0.282210 | 0.248877 | [0.1, 1.8] | No |
| First-birth fixed cost | 4.559138 | 4.559138 | [0, 8] | No |
| First-child room intercept | 0.364931 | 0.464931 | [0, 0.5] | No |
| 2007-2023 child-value change | -0.328714 | -0.328714 | [-1.5, 0.2] | No |
| Tenure dispersion $\kappa$ | 0.005000 | 0.005000 | Externally fixed | Not estimated |
| Dated supply elasticity $\eta$ | 0.630000 | 0.630000 | Externally fixed | Not estimated |

The old fertility-preference intercept is re-normalized to replacement for
every parameter vector. New candidate values are
$\psi_{2007}=0.28947879$ and
$\psi_{2023}=-0.03923512$.
The stored near-bound flag uses raw parameter units. For the bequest shifter,
interpret it alongside its position in the logarithmic search interval; a
raw-scale flag alone does not establish a binding search constraint.

The supply restriction 0.63 applies to dated markets. The inherited initialization
retains elasticity 1.75 before the dated supply intercept is normalized. Both
values and that order are unchanged across this comparison; whether this is the
intended economic contract remains a discussion item in the separate code audit.

The complete smaller-step panel has algebraic Jacobian rank
**11 of 11** and relative rank
**11 of 11** at the declared $10^{-3}$
threshold. These are local, step-dependent diagnostics, not standard errors or
proof of global identification. The saved receipt includes every singular value,
weak-direction loading, and forward/backward derivative comparison.
The weakest direction has **98.4%**
of its squared loading on the bequest shifter. The supply-intercept derivative
has forward/backward cosine **-0.369**;
opposite directional responses make a smooth local inverse unreliable here.

\clearpage

# Numerical checks and their limits

All completed cases pass the unchanged production acceptance checks. Those gates
do not certify every policy function. The separate adjacent-wealth screen flags
occupied value decreases in **2** evaluations; the largest exposed
pre-choice mass share is **2.08e-07**. These flags are retained and are
consistent with the already documented need to check the local saving optimizer;
this search does not establish their cause through a new global-optimizer solve.
The selected candidate has **0** flagged
occupied value decreases. Its reporting-floor budget-excess share is
**1.76e-09**. Probability arrays are finite
and remain inside $[0,1]$. The seventeen standard graphs follow this note unchanged.

The standard age-30 first-birth panel also shows a probability spike at 1 negative-wealth renter state. Its settled-childless pre-choice mass is 0. This is a check of current exposure, not a proof of optimal policies outside occupied states. The original graphs and their dense income-state labels are preserved.

The concurrent source-code audit also reproduces an incomplete policy-storage
interface: later-birth probabilities reside in a mutable parameter field instead
of the saved calendar policy. Reuse after an intervening solve can mix policies.
No affected result has been demonstrated in this historical calibration: the
nonzero preference drift clears warm starts at each later date. The interface
requires a separate repair and targeted replay before relying on vulnerable
policy-reuse paths. Exact reproduction alone would not detect this defect.

# What remains to discuss

The first-birth rooms contrast remains a provisional empirical-to-model mapping.
The young-ownership diagnostic is still excluded from estimation; model and
empirical age windows need a common contract, including existing prime-age and
old-wealth statistics. The overnight evidence on near-universal late-life
ownership and sensitivity to rental opportunities remains relevant. A better
in-sample objective does not independently resolve those issues.

Supply and credit results in the earlier morning audit belong to the retained
production calibration. This refinement supplies no revised policy magnitudes,
welfare ranking, or resolution of the household-entry transition.

The complete per-case fit and parameter tables are in each stage's
`report/all_target_fits.csv` and `report/all_parameters.csv` under
`output/model/e5f_morning_refinement_20260905a/`. Its README provides stage plans,
execution receipts, checkpoints, standard plots, and reproduction commands.
The selected case is **round2/task_009/summary.json**.

\clearpage

# Standard solution diagnostics

![Fertility by age.](output/model/e5f_morning_refinement_20260905a/round2/task_009/standard_diagnostics/fertility_by_age.png){width=88%}

![Ownership by age.](output/model/e5f_morning_refinement_20260905a/round2/task_009/standard_diagnostics/ownership_by_age.png){width=88%}

\clearpage

# Standard solution diagnostics

![Housing market.](output/model/e5f_morning_refinement_20260905a/round2/task_009/standard_diagnostics/housing_market.png){width=88%}

![Housing prices.](output/model/e5f_morning_refinement_20260905a/round2/task_009/standard_diagnostics/housing_prices.png){width=88%}

\clearpage

# Standard solution diagnostics

![Market clearing by market.](output/model/e5f_morning_refinement_20260905a/round2/task_009/standard_diagnostics/market_clearing_by_market.png){width=88%}

![Market clearing residuals.](output/model/e5f_morning_refinement_20260905a/round2/task_009/standard_diagnostics/market_clearing_residuals.png){width=88%}

\clearpage

# Standard solution diagnostics

![Owner rungs.](output/model/e5f_morning_refinement_20260905a/round2/task_009/standard_diagnostics/owner_rungs.png){width=88%}

![Tenure services.](output/model/e5f_morning_refinement_20260905a/round2/task_009/standard_diagnostics/tenure_services.png){width=88%}

\clearpage

# Standard solution diagnostics

![Fertility policy by age income state.](output/model/e5f_morning_refinement_20260905a/round2/task_009/standard_diagnostics/fertility_policy_by_age_income_state.png){width=88%}

![Housing by age income state.](output/model/e5f_morning_refinement_20260905a/round2/task_009/standard_diagnostics/housing_by_age_income_state.png){width=88%}

\clearpage

# Standard solution diagnostics

![Ownership by age income state.](output/model/e5f_morning_refinement_20260905a/round2/task_009/standard_diagnostics/ownership_by_age_income_state.png){width=88%}

![Liquid wealth by age income state.](output/model/e5f_morning_refinement_20260905a/round2/task_009/standard_diagnostics/liquid_wealth_by_age_income_state.png){width=88%}

\clearpage

# Standard solution diagnostics

![Income state outcomes.](output/model/e5f_morning_refinement_20260905a/round2/task_009/standard_diagnostics/income_state_outcomes.png){width=88%}

![Policy childless renter age30.](output/model/e5f_morning_refinement_20260905a/round2/task_009/standard_diagnostics/policy_childless_renter_age30.png){width=88%}

\clearpage

# Standard solution diagnostics

![Wealth dist childless renter age30.](output/model/e5f_morning_refinement_20260905a/round2/task_009/standard_diagnostics/wealth_dist_childless_renter_age30.png){width=88%}

![Policy childless renter age42.](output/model/e5f_morning_refinement_20260905a/round2/task_009/standard_diagnostics/policy_childless_renter_age42.png){width=88%}

\clearpage

# Standard solution diagnostics

![Wealth dist childless renter age42.](output/model/e5f_morning_refinement_20260905a/round2/task_009/standard_diagnostics/wealth_dist_childless_renter_age42.png){width=88%}
