---
title: "Overnight calibration: results and remaining gaps"
subtitle: "All eleven parameters searched under the unchanged twelve targets"
author: "Prepared for Tommaso De Santo"
date: "6 September 2026"
---

# The short read

**Best evaluated loss: 19.2844.** Final exact verification is incomplete; this is a provisional result.

The first-birth housing response is **0.5084 rooms against 0.7202**. Mean rooms are **6.3361 against 5.7800**, and ownership is **57.20% against 57.55%**.

| Point | Loss | First-birth rooms | Mean rooms |
|---|---:|---:|---:|
| Retained production | 30.4830 | 0.4364 | 6.3173 |
| Evening, old bounds | 21.7980 | 0.4766 | 6.3325 |
| Housing-only starting diagnostic | 20.6953 | 0.5158 | 6.3491 |
| Overnight joint search | 19.2844 | 0.5084 | 6.3361 |

All eleven estimated parameters could move during the overnight search. The first-child housing term could range from 0 to 2 instead of 0 to 0.5; every other bound, all twelve empirical targets and all weights stayed fixed. The model equations and numerical gates were retained.

The recovery stopped with status `failed`. There are 69 checked histories, including 4 smoke evaluations. Production remains the September 4 calibration; policy results have not been refreshed.

The first broad attempt stopped at candidate 17 when the housing-market residual exceeded the unchanged 0.0002 tolerance; none of its fourteen completed trials improved the seed. This separately recorded recovery uses smaller steps: each round probes all eleven parameters, then evaluates combinations of improving changes under the complete objective. It retains the original overnight cutoff.

**The economic result:** the joint search improves overall fit, chiefly on ownership and other margins. Its first-birth room response remains below the housing-only starting point. It has not closed the childbirth-housing gap.

**For the presentation:** retain the reproduced benchmark until the new point has passed exact repetition. The next numerical task is to diagnose the failed market solve before considering another search. The empirical target and all weights stay fixed.

\clearpage

# Complete target fit

Every gap is model minus target; the loss contribution is weight times squared gap. Provisional weight scales mean the total is an optimization objective, not a chi-squared statistic.

| Moment | Target | Model | Gap | Weight | Loss |
|---|---:|---:|---:|---:|---:|
| Completed fertility, age 42 | 1.9180 | 1.9276 | 0.0096 | 1425.74 | 0.1316 |
| Childless share, age 42 | 0.1880 | 0.1875 | -0.0005 | 17180.7 | 0.0039 |
| Age at first birth | 26.0446 | 26.3259 | 0.2813 | 44.4444 | 3.5160 |
| First births at age 30+ | 0.2603 | 0.2419 | -0.0184 | 10000 | 3.3807 |
| First-birth rooms response | 0.7202 | 0.5084 | -0.2118 | 137.565 | 6.1710 |
| Parent rooms: 3+ vs. 1-2 | 0.3677 | 0.3642 | -0.0035 | 2958.51 | 0.0364 |
| Parent ownership gap | 0.1677 | 0.1727 | 0.0050 | 14229.6 | 0.3577 |
| Ownership, ages 30-55 | 0.5755 | 0.5720 | -0.0035 | 1207.85 | 0.0147 |
| Mean rooms, ages 18-85 | 5.7800 | 6.3361 | 0.5562 | 11.9732 | 3.7034 |
| Wealth / annual earnings | 6.8731 | 7.0515 | 0.1784 | 6.28767 | 0.2002 |
| Bequests / wealth, annual | 0.0088 | 0.0085 | -0.0003 | 5.16529e+06 | 0.3695 |
| Old wealth/income: p90/p50 | 3.4481 | 3.2914 | -0.1567 | 56.9598 | 1.3995 |

Completed fertility and childlessness refer to the cohort centered at age 42 in 2023. The first-birth rooms response compares matched 2019-to2023 branches, allowing later births in the treated branch. Parent room gaps use ages 30--55. The final wealth statistic is p90/p50 of wealth divided by annual income at ages 76--84.

\clearpage

# Estimated parameters and restrictions

| Parameter | Estimate | Lower | Upper | Near bound? |
|---|---:|---:|---:|---|
| Annual patience | 0.995692 | 0.940000 | 0.999500 | No |
| First-birth dispersion | 2.238196 | 0.020000 | 50.000000 | No |
| Later-birth dispersion | 1.766446 | 0.020000 | 50.000000 | No |
| Owner premium | 1.049870 | 0.100000 | 5.000000 | No |
| Supply intercept | 14.427262 | 0.200000 | 80.000000 | No |
| Bequest strength | 0.566352 | 0.000000 | 8.000000 | No |
| Bequest shifter | 0.104376 | 0.020000 | 16.000000 | Yes |
| Per-child housing floor | 0.233741 | 0.100000 | 1.800000 | No |
| First-birth fixed cost | 4.555365 | 0.000000 | 8.000000 | No |
| First-child housing term | 0.606178 | 0.000000 | 2.000000 | No |
| Dated child-value change | -0.323676 | -1.500000 | 0.200000 | No |

Near-bound flags use the saved rule: within 2% of the raw parameter range. They need not imply proximity in the logarithmic search coordinate.

| Other recorded object | Value |
|---|---:|
| Old child-value intercept | 0.296278 |
| 2023 child-value intercept | -0.027397 |
| Tenure dispersion | 0.005000 |
| Dated supply elasticity | 0.630000 |

Dated supply elasticity 0.63 and tenure dispersion 0.005 are externally fixed. Old initialization retains elasticity 1.75; old fertility is normalized to 2.1. The dated child-value change is estimated. The inherited household-entry/population closure is retained and remains a limitation for interpreting long forecasts.

\clearpage

# Verification and evidence

**Why the recovery stopped.** At the first moving historical date (2011), round 2 joint proposal 7 failed the market check. The 18-step attempt ended at relative residual 0.0002596; the existing 30-step retry ended at 0.001174, above the unchanged 0.0002 acceptance limit. The failed proposal differs from the best completed point. No failed point was accepted, and no third search was launched.

The code had found a sign-changing price bracket and exhausted bisection without reaching the residual gate. The next diagnosis should inspect residuals and household choices inside that bracket. The saved log does not establish a discontinuity, coding error or nonexistence of equilibrium; it does not show that the empirical housing target is unreachable.

The source, target, history, market, accounting, feasibility and population checks passed for the selected point. Its occupied adjacent-wealth value screen flags **0 decreases**; the exposed pre-choice mass share is 0. Its reported-budget excess share is 1.605e-09.

The best point is not certified by two final exact repetitions. Inspect the failure/status receipt before using it.

All 69 completed cases passed a separate artifact and result audit: 1,449 original hashes checked on Torch. One case finished during cancellation, after the controller recorded 68; it is included in the final inventory and does not change the winner. The local collector checked 224 available original hashes.

Large checkpoints remain on Torch. The selected seventeen standard plots have been inspected. They retain irregular tenure policies and nearly universal ownership at old ages; a zero occupied-value-drop screen does not certify all household policies. Inherited diagnostic summary statistics use different clocks and age windows from the calibration table, which remains authoritative.

[Experiment design and reproducibility](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_full_joint_overnight_20260906a/recovery_01/README.md). [Completed inventory](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_full_joint_overnight_20260906a/recovery_01/completed_case_inventory.csv). [All target fits](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_full_joint_overnight_20260906a/recovery_01/final_target_fits.csv). [All parameter tables](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_full_joint_overnight_20260906a/recovery_01/final_parameters.csv).

The earlier local Jacobian had a weak bequest direction and a nonsmooth old wealth-to-income quantile. A lower loss alone does not resolve identification, age-window measurement, old-age ownership saturation or population-closure questions. These remain separate from the target-and-equation-preserving search.

\clearpage

# Selected standard diagnostics

These are four views from the unchanged seventeen-plot packet. Ownership below uses the full age distribution; the calibration row uses ages 30--55. The fertility panel is an inherited statistic of childless households' policies; the fitted completed-fertility target is a cohort stock.

![ownership by age](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_full_joint_overnight_20260906a/recovery_01/search/polish_2_joint/task_002/standard_diagnostics/ownership_by_age.png){width=48%} ![fertility by age](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_full_joint_overnight_20260906a/recovery_01/search/polish_2_joint/task_002/standard_diagnostics/fertility_by_age.png){width=48%}

![housing market](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_full_joint_overnight_20260906a/recovery_01/search/polish_2_joint/task_002/standard_diagnostics/housing_market.png){width=48%} ![housing prices](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_full_joint_overnight_20260906a/recovery_01/search/polish_2_joint/task_002/standard_diagnostics/housing_prices.png){width=48%}


The market plot checks the selected 2023 equilibrium, not the failed 2011 proposal. The old-age ownership concentration remains an economic-fit concern. Full policy, wealth-distribution, income-state and tenure panels accompany this report. Housing policy panels display a single selected branch; aggregate housing uses the tenure probabilities.

A visible age-42 fertility spike was traced to a wealth node with zero pre-choice and current childless-renter mass in the saved selected state. That checks its exposure at this solution; it is not a general guarantee about boundary policies.
