---
title: "Calibration: what improved and what needs a decision"
subtitle: "Evening follow-up for the September 14 presentation"
author: "Prepared for Tommaso De Santo"
date: "5 September 2026"
---

# The decision in one page

**There is a reproducible calibration improvement, but the main housing miss remains.** Within the existing parameter bounds, the new candidate lowers loss from 30.483 to **21.798**, a **28.5%** reduction versus retained production and **5.9%** versus the morning candidate. Two exact repeats confirm all twelve fit rows, parameters and the historical path.

Its first-birth housing response is **0.477 rooms against 0.720**: still a 33.8% shortfall. Mean housing remains **6.333 rooms against 5.780**, 9.6% too high. These discrepancies matter more than the improvement in the scalar objective.

| Point | Loss | First-birth rooms | Mean rooms |
| --- | ---: | ---: | ---: |
| Empirical target | -- | 0.7202 | 5.7800 |
| Retained production | 30.4830 | 0.4364 | 6.3173 |
| Morning candidate | 23.1534 | 0.4660 | 6.3286 |
| Evening, existing bounds | **21.7980** | **0.4766** | **6.3325** |
| Observed-gap diagnostic | **20.6953** | **0.5158** | **6.3491** |

The last row is a diagnostic fixed-parameter evaluation outside the old first-child bound; the remaining coordinates were not re-estimated. Two exact repeats confirm this point.


The model first-birth response compares rooms in 2023 between matched branches
originating with a first birth versus childlessness in 2019. Later births are
allowed in the birth branch; the control stays childless. The empirical target remains
**0.720246 rooms**. The earlier full-data regression check reproduced it and
gave no reason to lower it. Tonight's work keeps all twelve target values and
weights, the model equations and the production calibration intact.

**My recommendation:** keep the empirical target and pursue a joint housing-and-tenure refinement. The final diagnostic improves the first-birth response to **0.516** while matching the family-size room gap, at loss **20.695**. It provides a more promising direction than holding the three-child preference floor fixed. The mean-rooms and ownership gaps remain. Keep policy results tied to the morning calibration until a replacement is chosen and checked.

\clearpage

# What the housing experiment changes

The child housing floor is the minimum housing-service component in household
preferences. With $m>0$ children at home, it is

$$\bar h(m)=\bar h_J+m\bar h_C.$$

The first-child term $\bar h_J$ adds the same amount for every household with children at home;
the slope $\bar h_C$ adds more for each additional child. Raising the first term
by $d$ and lowering the slope by $d/3$ preserves the three-child floor and raises
the one-child floor by $2d/3$. It changes desired housing through an existing
model parameter.

The morning candidate has $\bar h_J=0.464931$ and $\bar h_C=0.248877$.
Its one-child floor is 0.713808 and its three-child floor is 1.211562.
The existing search ceiling $\bar h_J=0.5$ leaves little room for this particular
rebalancing. The ceiling applies to a preference parameter; equilibrium housing
responses also depend on saving, tenure, fertility choices and market prices.

The 0.5 ceiling entered the August 18 code when the contemporaneous estimate was about 0.20. The bounded history search recovered no empirical or theoretical justification for that endpoint and no direct author choice of it. The experiment README records the historical evidence.

The bounded search completed a 23-case local coordinate panel followed by twelve declared joint cases. It combines the first-child ceiling with alternative per-child floors and supply intercepts. The winning point preserves the three-child floor and combines half-steps along the improving non-housing coordinates.

A lower supply intercept can improve both mean rooms and the first-birth response, but reduces ownership. For example, one joint point gives a 0.502-room response and 6.059 mean rooms, with ownership 0.516 against the 0.575 target; its loss is 23.254. The current winner instead gives 0.477 rooms and ownership 0.547. These are observed trade-offs across the tested points, not a global frontier.

At the bounded winner, the matched birth branch uses **5.807702** rooms and the childless branch **5.331057**, giving the 0.476645 response. Relative to the morning point, the birth branch rises 0.007360 and the control falls 0.003242. Their difference raises the response by only **0.010601 rooms**. This comparison includes changes in equilibrium and the selected risk set; it is not a causal decomposition.

![The horizontal axis is the first-child term; the one-child floor equals this term plus the per-child slope. The three-child floor and nine other coordinates are fixed.](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_evening_housing_refinement_20260905a/analysis/conditional_floor_profile.png){width=100%}



\clearpage

# Complete fit under the retained search bounds

The best candidate within the existing bounds is the **evening round 2, case 10**, with loss **21.798033**. Its folder is distinct from the retained September 4 `task_010`. It uses the same twelve-row objective and eleven estimated parameters.

Every gap is model minus target. Loss contributions equal weight times squared
gap. Some weights use provisional scales, so the total loss is an objective
value, not a chi-squared test statistic.

| Moment | Target | Model | Gap | Weight | Loss |
| --- | ---: | ---: | ---: | ---: | ---: |
| Completed fertility, age 42 | 1.9180 | 1.9238 | 0.0058 | 1425.74 | 0.0485 |
| Childless share, age 42 | 0.1880 | 0.1883 | 0.0003 | 17180.7 | 0.0012 |
| Age at first birth | 26.0446 | 26.3433 | 0.2987 | 44.4444 | 3.9645 |
| First births at age 30+ | 0.2603 | 0.2429 | -0.0175 | 10000 | 3.0453 |
| First-birth rooms response | 0.7202 | 0.4766 | -0.2436 | 137.565 | 8.1633 |
| Parent rooms: 3+ vs. 1-2 | 0.3677 | 0.3727 | 0.0050 | 2958.51 | 0.0742 |
| Parent ownership gap | 0.1677 | 0.1679 | 0.0002 | 14229.6 | 0.0008 |
| Ownership, ages 30-55 | 0.5755 | 0.5465 | -0.0290 | 1207.85 | 1.0130 |
| Mean rooms, ages 18-85 | 5.7800 | 6.3325 | 0.5526 | 11.9732 | 3.6557 |
| Wealth / annual earnings | 6.8731 | 7.0374 | 0.1643 | 6.28767 | 0.1698 |
| Bequests / wealth, annual | 0.0088 | 0.0085 | -0.0003 | 5.16529e+06 | 0.4633 |
| Old wealth/income: p90/p50 | 3.4481 | 3.3031 | -0.1451 | 56.9598 | 1.1984 |
| **Total** | | | | | **21.7980** |

Shares are fractions; ages are years; the housing differences and mean rooms are in rooms. The parent room gap uses ages 30-55. The last row is the ratio of the 90th percentile to the median of wealth/annual income among ages 76-84. Weights are shown to six significant figures; the linked CSVs preserve full precision.

The first-birth housing response accounts for **37.4%** of this loss. Fertility timing contributes another **32.2%**, and mean rooms **16.8%**. Fertility levels and the parent ownership gap fit closely. This concentration suggests focusing the next search on housing demand and birth timing while retaining the full objective. The later-first-birth share fits better than at the morning point, while mean age fits worse; improving one timing statistic alone is insufficient.

\clearpage

# First diagnostic: the family-size trade-off

At fixed jumps 0.6, 0.75 and 0.9, first-birth room responses rise to 0.500, 0.540 and 0.580, while losses rise to **22.089, 30.679 and 49.496**. All three points pass the maintained gates and have no occupied value decreases. None improves the full objective. The anchor remains best and was already reproduced twice using the fixed option; two further identical repeats were unnecessary and were not launched.

Preserving the three-child preference floor still allows the observed large-family room gap to fall. At jump 0.9, this gap is **0.275 versus target 0.368**, adding 25.585 loss units relative to the anchor; the parent ownership gap adds another 5.641. The birth-rooms improvement does not offset these costs. This is evidence about the tested path, not proof that the ceiling is irrelevant or that joint re-estimation cannot help.

The complete fit below is for **jump 0.6**, the lowest-loss extended point. The other nine estimated coordinates stay fixed and the per-child slope compensates mechanically; this is not a ten-parameter recalibration.

| Moment | Target | Model | Gap | Weight | Loss |
| --- | ---: | ---: | ---: | ---: | ---: |
| Completed fertility, age 42 | 1.9180 | 1.9238 | 0.0058 | 1425.74 | 0.0486 |
| Childless share, age 42 | 0.1880 | 0.1889 | 0.0009 | 17180.7 | 0.0142 |
| Age at first birth | 26.0446 | 26.3476 | 0.3029 | 44.4444 | 4.0788 |
| First births at age 30+ | 0.2603 | 0.2431 | -0.0173 | 10000 | 2.9787 |
| First-birth rooms response | 0.7202 | 0.4995 | -0.2207 | 137.565 | 6.7011 |
| Parent rooms: 3+ vs. 1-2 | 0.3677 | 0.3477 | -0.0200 | 2958.51 | 1.1813 |
| Parent ownership gap | 0.1677 | 0.1727 | 0.0051 | 14229.6 | 0.3668 |
| Ownership, ages 30-55 | 0.5755 | 0.5457 | -0.0298 | 1207.85 | 1.0740 |
| Mean rooms, ages 18-85 | 5.7800 | 6.3425 | 0.5625 | 11.9732 | 3.7891 |
| Wealth / annual earnings | 6.8731 | 7.0354 | 0.1623 | 6.28767 | 0.1656 |
| Bequests / wealth, annual | 0.0088 | 0.0085 | -0.0003 | 5.16529e+06 | 0.4697 |
| Old wealth/income: p90/p50 | 3.4481 | 3.3017 | -0.1464 | 56.9598 | 1.2214 |
| **Total** | | | | | **22.0892** |


\clearpage

# Follow-up: fit the observed family-size room gap

The last four-point test uses the stable local housing-floor derivatives to choose a per-child slope that approximately preserves the **observed** parent-room contrast. It holds the other nine estimated coordinates fixed, re-normalizes old fertility, then evaluates the complete twelve-row objective. The best point is at a fixed first-child term **0.6** and slope **0.232688**.

It gives a first-birth response of **0.515778 rooms**, and a large-family room gap of **0.367541 against 0.367700**. Loss falls to **20.695274**, 5.1% below the bounded winner. This is a promising diagnostic point outside the original first-child bound; it is not a jointly re-estimated production calibration. **Both exact repeats passed.**

| Moment | Target | Model | Gap | Weight | Loss |
| --- | ---: | ---: | ---: | ---: | ---: |
| Completed fertility, age 42 | 1.9180 | 1.9238 | 0.0058 | 1425.74 | 0.0483 |
| Childless share, age 42 | 0.1880 | 0.1891 | 0.0011 | 17180.7 | 0.0208 |
| Age at first birth | 26.0446 | 26.3492 | 0.3046 | 44.4444 | 4.1230 |
| First births at age 30+ | 0.2603 | 0.2431 | -0.0172 | 10000 | 2.9525 |
| First-birth rooms response | 0.7202 | 0.5158 | -0.2045 | 137.565 | 5.7512 |
| Parent rooms: 3+ vs. 1-2 | 0.3677 | 0.3675 | -0.0002 | 2958.51 | 0.0001 |
| Parent ownership gap | 0.1677 | 0.1757 | 0.0081 | 14229.6 | 0.9250 |
| Ownership, ages 30-55 | 0.5755 | 0.5450 | -0.0305 | 1207.85 | 1.1240 |
| Mean rooms, ages 18-85 | 5.7800 | 6.3491 | 0.5691 | 11.9732 | 3.8784 |
| Wealth / annual earnings | 6.8731 | 7.0338 | 0.1607 | 6.28767 | 0.1624 |
| Bequests / wealth, annual | 0.0088 | 0.0085 | -0.0003 | 5.16529e+06 | 0.4735 |
| Old wealth/income: p90/p50 | 3.4481 | 3.3008 | -0.1473 | 56.9598 | 1.2361 |
| **Total** | | | | | **20.6953** |

The first-child term is fixed above its original [0, 0.5] interval. The per-child slope is inside [0.1, 1.8] and not near either bound. All other nine coordinates and their restrictions are exactly those in the next page's bounded-parameter table. The normalized old child-value intercept is **0.299775** and its derived 2023 value **-0.027106**. No target or weight changed.

Compared with the bounded point, the birth branch gains 0.0465 rooms and the childless branch gains 0.0073, raising their difference by 0.0391 rooms. This comparison includes equilibrium and risk-set changes.

The gain has limits: mean rooms rise to 6.349, and the parent ownership gap widens to 0.1757 against 0.1677. The first-birth response still falls short by **0.2045 rooms**. This supports joint housing-and-tenure refinement while retaining all twelve moments; it does not close the calibration problem.

\clearpage

# All estimated parameters and restrictions

| Parameter | Bounded estimate | Search interval | Near bound? |
| --- | ---: | --- | --- |
| Annual patience $\beta$ | 0.995739 | [0.94, 0.9995] | No |
| First-birth dispersion $\kappa_1$ | 2.232730 | [0.02, 50] | No |
| Later-birth dispersion $\kappa_C$ | 1.736471 | [0.02, 50] | No |
| Owner premium $\chi$ | 1.043472 | [0.1, 5] | No |
| Supply intercept $H_0$ | 14.562959 | [0.2, 80] | No |
| Bequest strength $\theta_0$ | 0.554442 | [0, 8] | No |
| Bequest shifter $\theta_1$ | 0.106357 | [0.02, 16] | Raw scale |
| Per-child floor $\bar h_C$ | 0.237187 | [0.1, 1.8] | No |
| First-birth fixed cost | 4.544053 | [0, 8] | No |
| First-child term $\bar h_J$ | 0.500000 | [0, 0.5] | Upper ceiling |
| Dated child-value change | -0.326881 | [-1.5, 0.2] | No |

In the three-child-floor-preserving tests, **all nine coordinates other than the two housing floors retain the values above**. The first-child term is externally fixed at the diagnostic value; the slope is set to preserve the three-child floor. Neither is estimated in this exercise. The slope at jump 0.9 is near its original lower bound; the first-child values are outside the original [0, 0.5] interval.

| Fixed first-child term | Per-child slope | Normalized old intercept | Derived 2023 intercept |
| ---: | ---: | ---: | ---: |
| 0.600000 | 0.203854 | 0.297403 | -0.029477 |
| 0.750000 | 0.153854 | 0.301649 | -0.025232 |
| 0.900000 | 0.103854 | 0.305974 | -0.020907 |


The near-bound flag follows the stored rule: within 2% of the range in raw
parameter units. The bequest shifter can therefore be flagged even when it is
well inside its logarithmic search coordinate. The dated housing-supply
elasticity is fixed at 0.63 and tenure-choice dispersion at 0.005. The inherited
old-equilibrium initialization still uses elasticity 1.75. The old fertility
intercept is normalized to completed fertility 2.1; the dated preference change
is an estimated coordinate. None of these restrictions is newly estimated here.

At the bounded winner, the normalized old child-value intercept is 0.294616 and the derived 2023 intercept is -0.032265. The young-ownership validation row remains outside the objective: the retained whole-node convention gives 31.094% versus the ACS 34.117% target. The coarse age-cell measurement is recorded separately in the diagnostic packet.

\clearpage

# What the numerical checks establish

The evening work completed **50 full historical evaluations**: 39 in the bounded search, five in the first diagnostic and six in the observed-room-gap follow-up. The same source, target fingerprint, population/accounting contract, market tolerance and feasibility checks pass. Each history has a 2023 checkpoint and the unchanged seventeen-graph packet. There were also two fixed-price Bellman evaluations for the value warning described below. Allocated cluster time was **8.05 CPU-hours**.

The collector verified all **1,050 originally recorded artifact hashes** on Torch. All **2,742 collected local files** match the final manifest; **52 large checkpoints** remain on Torch. The bounded final repeats, fixed-option smokes and final diagnostic repeats reproduce all twelve fit rows, all physical parameters, 253 numeric historical entries and all seventeen standard PNGs exactly. The five launcher/planning guard checks pass.

The bounded winner and all seven extended diagnostic points have zero occupied value decreases. Maximum reported-budget excess share across the histories is $1.19\times10^{-8}$. Inspection of both selected points and the largest extension preserves the familiar equilibrium shapes: ownership still approaches one at old ages, and tenure policies remain uneven near some wealth thresholds. A lower objective does not settle these validation issues.

The coordinate panel around the morning starting point gives a weighted
Jacobian with condition number **1,335**. It has eleven nonzero singular values,
but only ten exceed the relative $10^{-3}$ screening threshold. Its weakest
direction is dominated by the bequest shifter. This is a warning about local
conditioning under the chosen step and parameter scaling; it is not proof of
global underidentification or an identification certificate for the new winner.

The main disagreement between upward and downward derivatives comes from the
old-age wealth-to-income p90/p50 target. Nine saved-state reconstructions reproduce that
statistic exactly. The code uses an inverse-CDF quantile on discrete support:
a small discount-factor reduction moves p90 from 27.895 to 27.587 while p50
stays near 8.443. This explains an objective jump without implying a coding
error. The housing-floor derivatives are much more consistent. No moment was
dropped, smoothed or reweighted. The full Jacobian was not inverted; the final four-point diagnostic used only the stable floor derivatives to choose explicit test points.

One rejected tenure-preference case has a value decrease over an adjacent
wealth step, affecting 0.229 per million of pre-choice mass. A fresh local
solve reproduces its policies and distributions exactly. Exhaustive saving
removes the decrease; at the same prices and distribution, births change by
0.000234% and ownership by -0.000022 percentage points. This does not bound
the error of a re-solved historical calibration or of the matched birth-rooms
statistic. The verified quantity differences are too small to account for the
large ownership movement at that saved state.

\clearpage

# Policy context and next discussion

These policy numbers were verified at the **morning candidate, loss 23.153**.
They have not been recomputed at tonight's new points. Entries are percent
changes in births per household against that candidate's own no-policy path.

| Policy | Impact, 2023 | 2063 |
| --- | ---: | ---: |
| Supply +20% | +1.343% | +2.141% |
| Dependent-child LTV 95% | +0.018% | +0.075% |
| Supply and credit | +1.366% | +2.244% |
| Property tax 1% to 2%, no rebate | -1.145% | -1.699% |

The separate rebated 2%-versus-rebated 1% impact decomposition is -1.377%
direct tax, +0.302% capitalization and +1.574% rebate, totaling +0.500%.
Its fiscal baseline differs from the unrebated tax row. These remain conditional
mechanism comparisons, with an inherited household-entry closure; they are
neither welfare results nor validated population forecasts.

The next discussion should settle three things: whether to enlarge the first-child parameter's search interval on the evidence above; which candidate is ready for the September 14 presentation; and which population-closure assumptions can support the policy interpretation. Keep the full twelve-row objective in the next search. The current evidence does not establish that the first-birth response is unreachable.

The empirical reference-choice issue remains a documentation item. The completed paired regression reproduces 0.720246 rooms; changing one cohort reference gives 0.730459 with identical observations and fitted values. This demonstrates sensitivity of the aggregation convention, but the tested change is small and upward. It supplies no basis for lowering the calibration target. The pooled empirical cohort mix and matched model branch should also be described distinctly.

# Evidence and reproduction

The [experiment README](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_evening_housing_refinement_20260905a/README.md) records the design, deadlines, jobs, immutable source and target hashes, and graph replay command. The [complete bounded-case ledger](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_evening_housing_refinement_20260905a/analysis/all_cases.csv), [all twelve-row fits](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_evening_housing_refinement_20260905a/analysis/all_target_fits.csv) and [all parameter tables](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_evening_housing_refinement_20260905a/analysis/all_parameters.csv) include the rejected points.

Earlier comparison points retain their complete tables: [production fit](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_transition_calibration_eta063_kappa005_rooms_repair_gauss_newton_coord_20260904a/task_010/target_fit_long.csv), [production parameters](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_transition_calibration_eta063_kappa005_rooms_repair_gauss_newton_coord_20260904a/task_010/parameter_table.csv), [morning fit](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_morning_refinement_20260905a/round2/task_009/target_fit_long.csv), and [morning parameters](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_morning_refinement_20260905a/round2/task_009/parameter_table.csv). The [earlier policy review](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/docs/model/e5f_candidate_policy_comparison_review.md) contains the full policy evidence. The [empirical measurement review](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/docs/model/e5f_first_birth_measurement_review.md) records the regression replay and its limits.

The [conditional profile ledger](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_evening_housing_refinement_20260905a/analysis/profile_cases.csv), [complete profile fits](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_evening_housing_refinement_20260905a/profile_grid/report/all_target_fits.csv), [profile parameters](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_evening_housing_refinement_20260905a/profile_grid/report/all_parameters.csv) and [stop receipt](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_evening_housing_refinement_20260905a/profile_stop.json) record the extended points and why further identical repeats were omitted.

The [final four-point ledger](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_evening_housing_refinement_20260905a/empirical_room_gap/grid/report/all_candidates.csv), [full fits](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_evening_housing_refinement_20260905a/empirical_room_gap/grid/report/all_target_fits.csv) and [parameter ledger](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_evening_housing_refinement_20260905a/empirical_room_gap/grid/report/all_parameters.csv) record the observed-room-gap follow-up. Its preparation and corrected collector are saved in the experiment folder.

The standard seventeen diagnostic plots accompany every completed historical
evaluation. The saved-state graph replay reproduced all seventeen original PNG
hashes exactly with no model solve. The experiment README gives the single
regeneration command, all plans, budgets and source restrictions. Production
remains the September 4 `task_010`; author-controlled manuscript files were not
edited.
