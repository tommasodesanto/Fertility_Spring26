---
title: "First-birth housing: a measurement issue to resolve"
subtitle: "Calibration assessment for the September 14 presentation"
author: "Prepared for Tommaso De Santo"
date: "5 September 2026"
---

# The decision for discussion

**Resolve the empirical housing statistic before another calibration search.**
The audit identifies a normalization problem in the event-study design: some
birth cohorts have no observation in the stated reference period, and their
weights differ between the two endpoints used to construct the target. An
arbitrary shift in their coefficient levels can therefore change the reported
contrast without changing fitted housing outcomes. An independent read-only
review confirms the algebra. The exact omitted coefficients in the production
regression still need to be recovered.

This finding does **not** supply a corrected rooms target, establish the size or
sign of its error, or explain the whole model shortfall. The recorded target
remains **0.720246 rooms**, with its existing reported standard error **0.085260**.
The existing model, target file, weights and calibration results are unchanged.
The appendix retains every target-fit row and estimated parameter for both the
production calibration and the verified refinement candidate.

# Why the reference period matters

A birth cohort is the group of women with a first birth in the same calendar
year. The estimator fits a separate housing profile for each cohort and then
averages those profiles using the cohort shares observed at each event time.
The reference period, two years before first birth, is assigned a zero effect.

The sample-only replay finds **16 birth cohorts with no observation at event
year -2**. Such cohorts account for **41.91%** of the cohort-share weight at -1
and **46.01%** at +3. A second sample-only check verifies that all regression
inputs are nonmissing, first-birth cohort is constant within person, and the
included event indicators cover every treated observation except -2. These
checks require no regression coefficients.

For a cohort without the reference observation, shifting all its observed
housing coefficients up by one room and its person fixed effects down by one
room leaves every fitted observation unchanged. If that cohort receives
different weights at the two endpoints, the weighted difference changes.
That is a failure of invariance to an arbitrary normalization, beyond the usual
question of whether changing cohort composition is economically desirable.

\clearpage

# A concrete check of the normalization

Let $\delta_{g,k}$ denote the cohort-$g$ housing coefficient at event time $k$,
and $w_{g,k}$ its cohort share at that event time. The reported target is

$$
D=\sum_g w_{g,3}\delta_{g,3}-\sum_g w_{g,-1}\delta_{g,-1}.
$$

For a cohort with no -2 observation, the transformation
$\delta_{g,k}\mapsto\delta_{g,k}+c_g$ on its supported event cells and
$\alpha_i\mapsto\alpha_i-c_g$ for its members' person fixed effects leaves the
regression's predictions unchanged. The reported contrast, however, changes by

$$
\Delta D=(w_{g,3}-w_{g,-1})c_g.
$$

The **1986 first-birth cohort** supplies a numerical example with observations
at both endpoints. It has no -2 observation; its weights are **0.026299** at -1
and **0.040435** at +3. A one-room normalization shift changes the contrast by
**0.014136 rooms**, while the verified change in every fitted row is **zero**.
The one-room shift is an algebraic illustration, not an economic intervention,
an estimate of bias, or a proposed correction.

The installed estimator calculates these cohort shares before the regression's
singleton pruning. The all-inputs-observed assertion verifies that the saved
input support is the sample used for that share calculation. The zero-prediction
identity also holds on any subset retained by the regression. What remains
unrecovered is the realized estimation sample, the omitted-column choices and
the cohort coefficient/covariance objects needed to quantify an alternative
well-defined contrast. The original receipt reports 49,457 estimation rows;
the verified input sample contains 49,872 rows.

A common set of cohort weights applied to estimable within-cohort changes would
cancel this particular arbitrary level. That is a candidate diagnostic repair,
not an approved replacement target: cohort overlap, remaining rank conditions,
uncertainty and economic interpretation must all be checked. No target is dropped,
demoted or reweighted in this review.

# What passed, and what did not

The source and target provenance verification passed, including the raw PSID
source hash. The first sample-only replay passed in **28 seconds**, and the
additional design assertions passed in **84 seconds**. The one predeclared full
regression replay reached its **900-second limit** and was stopped; its recorded
elapsed time is 902.5 seconds including termination. It produced no usable
cohort coefficients or covariance matrix. No second regression was launched.
Post-run hashes confirm that the primary builder and every primary target
output remain unchanged. Only cohort/event aggregates were exported.

\clearpage

# The event-time pattern does not supply a convenient replacement

![Supplemental empirical profile from the unchanged production event table. Bars show its existing pointwise 95 percent intervals, conditional on the reported specification and normalization.](output/model/e5f_first_birth_measurement_review_20260905a/empirical_event_profile.png){width=100%}

The profile alternates substantially between odd and even event years. However,
nearby four-year summaries of the saved coefficients remain near 0.72 rooms:

| Contrast | Rooms | Interpretation |
|---|---:|---|
| -1 to +3 | 0.720246 | Existing target; reported SE 0.085260 |
| -2 to +2 | 0.726160 | Diagnostic; reported SE 0.080315 |
| Mean of +2/+3 minus mean of -2/-1 | 0.723203 | Diagnostic; joint uncertainty unavailable |
| 0 to +4 | 0.385433 | Different object, excluding adjustment through the birth year |

These comparisons are arithmetic diagnostics from the old coefficient table,
not newly estimated or validated causal effects. The missing covariance prevents
an uncertainty calculation for their differences. In particular, the smaller
0-to-4 value is not a reason to replace the target with a number closer to the
model.

The treated input observations span **1970-2019**. The weighted mean first-birth
year is **1994.64** at -1 and **1992.97** at +3; the share of weight on birth
cohorts from 2007 onward falls from **21.39% to 18.05%**. The total variation
distance between endpoint cohort shares is **0.20236**: 20.24 percentage points
of cohort probability weight would have to be reassigned to make the two mixes
identical. These are cohort-share input statistics, not claims about the final
regression sample or a quantified correction to the target.

\clearpage

# What the model measures

The active housing statistic compares two copies of the households that successfully
have a first birth. One copy has the baby; the other remains childless. Both begin
with the same wealth, age, income, tenure and location, and follow their respective
saving and housing choices over one four-year model period. The difference in their
destination rooms is also a difference in changes because their starting states
are identical.

At the review candidate, the 2019-2023 treated copy occupies **5.800342 rooms**
and the childless copy **5.334298**, giving **0.466044 rooms**. The treated branch
already permits a later birth at the destination date: continuation births equal
**26.88%** of its mass. Therefore, a missing second-birth option in this active
measurement is not the explanation for the shortfall. An older same-policy helper
has different timing and must not be confused with this dated calculation.

| Comparison | Empirical target | Model statistic |
|---|---|---|
| Population | PSID women who are reference persons or spouses, with one observation per household-year | Successful first-birth households in the fitted first-birth risk set |
| Control | Women confirmed childless in their observed histories | Identical copies held childless through the window |
| Clock | First-birth event years -1 and +3 | One four-year interval, 2019-2023 for the fitted terminal moment |
| Calendar coverage | Pooled historical PSID interviews | The fitted terminal four-year cohort |
| Before-birth differences | Estimated from an observational panel; anticipation and selection remain relevant | Starting-state differences between copies are zero by construction |
| Later births | Allowed in treated mothers' histories | Allowed in the treated branch at the destination |

Matching the duration does not establish that these are identical economic
objects. In particular, the model comparison holds the initial state fixed for
households selected into first birth, whereas the empirical estimator uses a
separate childless population. Neither exact reproduction nor a small standard
error establishes the missing comparability assumptions. No new estimate of the
size of this model-to-data selection discrepancy is available from this replay.

\clearpage

# What this means for calibration

The best review candidate reduces the unchanged twelve-target loss from
**30.482967 to 23.153445**. Across the **35 coordinate and joint evaluations**,
the largest first-birth response is **0.466044 rooms** and the smallest average
housing level is **6.270186 rooms**. Those extrema do not occur at the same
parameter vector. The target pair is **0.720246** and **5.779970 rooms**.
This neighborhood contains a meaningful improvement but provides no evidence
that the two housing discrepancies have been jointly resolved. It is not a
proof that the target pair is globally unreachable.

![Supplemental housing-fit comparison using saved evaluations. No new model solution is computed.](output/model/e5f_first_birth_measurement_review_20260905a/observed_housing_fit.png){width=100%}

After resolving the empirical normalization, the model measurement exercise should reproduce a stated empirical
population and event clock in the model, including before-birth housing choices,
and compare that statistic with the current matched-copy response. It should
preserve the full twelve-moment fit table and distinguish a measurement change
from an improvement under the existing objective. A new empirical weighting
rule requires its own uncertainty calculation and target-contract decision.

Until that comparison is resolved, retain both the production calibration and
the verified refinement candidate for discussion. The principal question is how
much of the rooms shortfall reflects housing-demand parameters, and how much
reflects comparing different populations and timing conventions. Another broad
search alone cannot answer it. Existing policy numbers belong to the retained
production calibration; this report supplies no revised policy or welfare result.

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

# Evidence and reproduction

All paths below are relative to the project root:
`/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26`.

- Primary estimator: `code/data/psid_followup_mar2026/sa_rooms_first_birth_household_aligned_v1.do`, sample construction and event definitions at lines 65-147, estimator and contrast at lines 166-212.
- Primary target and original event table: `code/data/psid_followup_mar2026/output/sa_rooms_first_birth_household_aligned_v1/`.
- Sample replay driver: `code/data/psid_followup_mar2026/audit_first_birth_rooms_measurement.py`.
- Aggregate analysis and two supplemental figures: `code/model/tools/analyze_e5f_first_birth_measurement.py`.
- Complete receipts, aggregates, null-direction check, failed regression receipt and independent review: `output/model/e5f_first_birth_measurement_review_20260905a/`.
- Active matched-branch implementation: `code/model/tools/run_e5f_transition_calibration.py`, functions `begin_dated_first_birth_housing_branch` and `finish_dated_first_birth_housing_branch`.
- Saved candidate: `output/model/e5f_morning_refinement_20260905a/round2/task_009/`; its dated-branch CSV is under `cases/task_001/` within that folder. This inner case identifier does not change the selected outer candidate.
- Earlier complete calibration readout and standard seventeen-graph packet: `docs/model/e5f_bounded_calibration_refinement_review.md` and `output/pdf/e5f_bounded_calibration_refinement_review.pdf`.

The installed estimator was also inspected at
`/Users/tommasodesanto/Library/Application Support/Stata/ado/plus/e/eventstudyinteract.ado`:
cohort-share estimation at lines 40-47, interacted regression at line 85, and
weighted aggregation at lines 103-108. The evidence manifest pins its hash.

Regenerate the aggregate analysis without a Stata regression or model solve:

```sh
python3 code/model/tools/analyze_e5f_first_birth_measurement.py
```

The output-folder README retains the original regression budget, both successful
sample checks, the failed full replay and the PDF build command. Full per-case
target and parameter tables are retained alongside this report. The two new
figures are supplemental; the standard model diagnostic packet is unchanged.
