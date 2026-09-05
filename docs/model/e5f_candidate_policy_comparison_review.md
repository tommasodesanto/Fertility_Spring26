---
title: "Calibration and policy: a clear comparison"
subtitle: "Decision readout for the September 14 presentation"
author: "Prepared for Tommaso De Santo"
date: "5 September 2026"
---

# The result to take into our discussion

**The verified candidate improves calibration fit and preserves the policy ordering.**
The unchanged twelve-target loss falls from **30.483 to 23.153**,
a **24.04%** improvement. The largest change in
the reported impact or 2063 birth response is **0.090 percentage points**.
The candidate remains a review point; the retained production calibration has
not been overwritten.

Supply expansion increases housing use and births. Family mortgage credit
increases ownership with little change in housing use or births. An unrebated
property-tax increase reduces births. Equal tax rebates alter that conclusion
because they change household resources. The fiscal convention must accompany
every tax claim.

**The main calibration weakness remains visible:** first-birth rooms are
**0.466 versus 0.720** in the data; mean rooms are **6.329 versus 5.780**.
The completed empirical check supports retaining the 0.720 target. A better
objective value does not resolve these housing discrepancies or establish that
the policy fertility responses are directly identified by causal cost variation.

My recommendation is to use this candidate as the working presentation point
after our review, with these limitations explicit. Keep the current benchmark
as the documented comparison. Another broad calibration search is not the
immediate priority.

# How much the policy numbers change

Entries below are percent changes in **births per household** against each
calibration's own no-policy path. 2063 is a finite transition date.

| Policy | Current 2023 | Candidate 2023 | Current 2063 | Candidate 2063 |
| --- | --- | --- | --- | --- |
| Supply +20% | +1.291 | +1.343 | +2.050 | +2.141 |
| Family LTV 95% | +0.018 | +0.018 | +0.076 | +0.075 |
| Supply and family credit | +1.314 | +1.366 | +2.155 | +2.244 |
| Tax 1% to 2%, no rebate | -1.100 | -1.145 | -1.623 | -1.699 |

\clearpage

# What moves on impact

| Candidate policy | Housing user cost | Rooms / HH | Ownership | Births / HH |
| --- | --- | --- | --- | --- |
| Supply +20% | -18.192% | +5.743% | +1.836 pp | +1.343% |
| Family LTV 95% | +0.060% | +0.044% | +0.732 pp | +0.018% |
| Supply and family credit | -18.017% | +5.888% | +2.988 pp | +1.366% |
| Tax 1% to 2%, no rebate | +14.577% | -4.923% | -2.334 pp | -1.145% |

The supply experiment shifts the supply schedule by 20% at every given price.
Equilibrium prices adjust, so the resulting increase in occupied rooms is the
smaller amount reported above. Family credit applies to households with
dependent children, not only first-time buyers.

For a concrete scale, 10,000 model households have about **846.1 births**
in the baseline four-year period. Supply expansion adds about
**11.4 births**; family credit adds about
**0.15**. This is an illustration of model units,
not an annual national forecast.

![Supplemental comparison of the complete policy paths. Dashed grey is the retained benchmark; solid blue is the candidate. The standard seventeen model diagnostics are preserved separately for every candidate policy date.](output/model/e5f_candidate_policy_comparison_20260905a/policy_fertility_comparison.png){width=100%}

\clearpage

# Tax results depend on the rebate comparison

The following three-factor Shapley decomposition averages each channel's
contribution across all six orders in which tax, house price and rebate can
change. Entries are percent of births per household in the **rebated 1%**
baseline; the reform is **rebated 2%**. Both equilibria use the same inherited
2023 distribution. The eight component cells reproduce the endpoints and add up.

| Birth contribution | Current benchmark | Verified candidate |
| --- | --- | --- |
| Direct tax | -1.323% | -1.377% |
| House-price capitalization | +0.291% | +0.302% |
| Equal rebate | +1.531% | +1.574% |
| Net | +0.498% | +0.500% |

This net effect is not the effect relative to the unrebated status quo. The
five-policy table instead compares unrebated 2% with unrebated 1%. The earlier
September 4 chat explicitly drew this distinction; the present review preserves it.

# Birth rates and total births in 2063

| Candidate policy | Births / HH | Household mass | Total births |
| --- | --- | --- | --- |
| Supply +20% | +2.141% | +0.418% | +2.568% |
| Family LTV 95% | +0.075% | +0.008% | +0.082% |
| Supply and family credit | +2.244% | +0.428% | +2.682% |
| Tax 1% to 2%, no rebate | -1.699% | -0.327% | -2.020% |

Total births combine the birth rate and the number of households. The exact
identity is checked for every row. Household mass is identical across policies
before 2043, as required by the maintained twenty-year entry lag. The inherited
household queue is a conditional mechanism closure; this table is not a
validated resident-population forecast or a terminal steady-state comparison.

\clearpage

# What was checked, and what remains open

All five short transition tests passed before the full five paths were launched.
Every full path exactly reproduces its short test's first two dates. The historical
reconstruction, shared starting state, market clearing, population, feasibility,
probability and accounting checks pass at the existing tolerances. Each of the
55 policy dates has a saved state and the unchanged seventeen-graph packet.

The source contains the already verified policy-reuse correction. Its fixed-
parameter history reproduces all twelve target rows, every parameter and bound,
and 253 numeric historical entries of the original candidate exactly. No model
primitive, target, weight, parameter or numerical tolerance was changed here.

The numerical screen asks whether household value ever falls when wealth rises.
It flags **28 of 55 dates**, but the largest affected share is only
**1.77 per million** of the inherited household distribution.
For three flagged policy-date distributions, a six-solve check first reproduces
the saved policies exactly, then searches the full saving grid at the same prices.
This removes the flagged value decreases in all three tests. The largest absolute
change in births per household is **0.000346%**. The eleven tax
equilibria/component cells contain 4 flagged packets, with maximum
affected share **2.58 per million**. These checks do not bound
the error of a fully re-solved equilibrium path or the tax decomposition.
The largest share affected by the reporting consumption floor exceeding the
budget is **0.022 per million**. Both rebated equilibria also pass
the check that the selected rebate matches the parameter state used for accounting.

The old chats explain several apparent inconsistencies. Supply elasticity 1.75
initializes the old equilibrium; 0.63 applies to dated markets. This distinction
was already recorded, though its economic justification remains open. A separate,
later demographic branch tracks people more consistently; it does not change
this retained household-entry experiment. The rebated-tax result uses its own
fiscal baseline. Direct author statements support studying policy during the
transition, without settling every inherited population assumption.

For the presentation, the defensible claim is that these mechanisms survive a
meaningful improvement in fit under the same target system. The remaining work
is to decide how prominently to qualify the housing-fit shortfall, population
closure and fiscal scope. This readout does not certify a new perfect-foresight
equilibrium, welfare ranking or globally optimal calibration.

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


\clearpage


# Reproduction and evidence

The evidence folder is `output/model/e5f_candidate_policy_comparison_20260905a`. It contains the pinned plan,
historical reconciliation, source and input hashes, all policy rows, exact
smoke/full comparison receipts, tax decomposition, numerical-flag check and
complete fit/parameter tables. Its README states the run budgets and commands.
The new comparison CSV contains all 110 policy/year/calibration rows, including
the no-policy paths. The standard graphs are in `full_diagnostics/`, grouped by
policy and date. All paths are relative to:
`/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26`.

Rebuild this readout from the collected outputs:

```sh
python3 code/model/tools/build_e5f_candidate_policy_comparison_report.py
```
