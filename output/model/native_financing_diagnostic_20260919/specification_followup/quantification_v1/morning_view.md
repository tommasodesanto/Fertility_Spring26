# Quantification: proposed baseline and the next calibration

**Working lead recommendation, September 20 evening.** Numerical-panel results
are pending. This document must be updated from actual smoke/production receipts
before the morning report. It is a proposal for the author, not an adopted
baseline or a claim that the calibration is complete.

The next milestone is a coherent, inspectable baseline. The hard tenure/space
restriction is a legitimate approximation to limited access to family-sized
rentals. Its literal support need not match every observed dwelling to be useful.
The relevant test is whether a parsimonious housing block can jointly account
for space, ownership, family differences and lifecycle profiles with credible
parameters. Policy responses from the current underfitting checkpoints cannot
settle that question.

## Recommended working choices

| Object | Recommendation | Quantification task before a serious refit |
|---|---|---|
| Earnings architecture | Use an age profile, persistent AR(1) and iid transitory shock, without a permanent type in the proposed baseline. Keep a fixed-effect process as a named robustness alternative. | Document the PSID gross head-plus-spouse earnings definition and age/sample restrictions; report the no-type long-lag covariance misfit. Sharing a functional form with BGM does not justify copying its income variances. |
| Income concept and taxes | Keep gross labor earnings with the model's explicit payroll tax for the working proposal. | Align the earnings process, age profile and wealth/earnings denominator. Treat measurement error through a sourced assumption or an explicit sensitivity, never an arbitrary variance factor chosen for model fit. |
| Four-year income and state grid | Prefer a moment-matched period-average proxy as the candidate numerical approximation; do not freeze 15 states yet. | Validate period covariance approximation, discrete level dispersion and conditional choices. Separate the change in discretization from changes in entrant wealth-income composition. The existing 15/27/45 evidence is a joint numerical/entry sensitivity. |
| Housing access and ownership preference | Keep the hard rental-size restriction and uniform ownership-service preference as the reference. | Fit housing levels and parent/nonparent differences jointly. The positive rental-cost tests did not identify a superior replacement. A family-specific ownership preference needs separately informative evidence before becoming another free parameter. |
| Fertility choices | Keep the retained sequential attempt/realization/tenure architecture, fecundity schedule and existing fertility taste scales. | Inspect how the existing coordinates jointly determine timing, childlessness and housing responses. Do not change information timing or add preference heterogeneity to obtain a desired policy sign. |
| Child space, goods costs and children at home | Keep the existing space requirement, cost terms and resident-child approximation provisionally. | Compare children at home and housing by parent age under matched definitions before changing departure hazards or adding an earnings penalty. Decreasing marginal components alone are not evidence that total child costs are wrong. |
| Mortgage, bequests, tenure smoothing | Keep the existing contract, estate valuation and fixed tenure scale for the controlled calibration diagnosis. | Record these as explicit maintained restrictions. Full amortization needs a defensible asset/principal representation; estate liquidation needs a receiver/valuation contract; deterministic tenure needs a market-clearing rule. These are separate extensions. |
| Target geography | Keep the existing 42-metro housing contract for all overnight comparisons. | Before adopting a new baseline, state whether the quantitative object is a metropolitan economy with national external inputs or a nationally matched population. Do not switch geography because its target values are easier to fit. |

This freezes enough structure for useful diagnosis while separating empirical
inputs from free structural parameters and numerical choices. It does not claim
that all those provisional restrictions have been validated.

## What the existing search can establish

The saved 96 proposals contain 89 complete valid fits and seven rejected
proposals. The lead independently reproduced every valid weighted loss from
the 12 scored rows and checked the complete target/weight signature. The
separate thirteenth row is the maintained stationary fertility normalization.
Full tables, actual bounds and descriptive trade-offs are in
[saved_fit/](saved_fit/README.md); the previously selected point's full tables
and exact repetitions are in
[the original readout](../../overnight/final_search/readout.md).

These are observations from a finite adaptive search. They do not establish an
attainable fit frontier or show that the current architecture cannot fit the
housing and fertility rows. We therefore did not adopt Claude's proposed family
ownership premium or its arbitrary fit cutoffs. See the
[lead review](claude_review/lead_review.md).

The four evaluated points closest jointly to the mean-room and overall
ownership targets have poor overall fit, chiefly because their recent-parent
ownership gap is about 0.23–0.25 below its target. This is a joint-fit problem
in the observed sample, not evidence that either level is individually
unreachable. Full rows remain in the linked tables. The large contribution
reflects the frozen working weight as well as the gap.

That row deserves explicit measurement attention before a new preference is
added. The active scorer uses the approved synchronized proxy: current births
from previously empty-dependent homes versus current empty homes, including
former parents, ages 30–55. ACS instead compares households whose oldest
resident own child is under four with households having no resident own
children. This is an explicitly maintained approximation, not an unnoticed
substitution of all parents for recent parents. The stored older warning about
a lifetime-childless control does not describe the active synchronized observer.
Sampling uncertainty alone does not measure approximation error. Any later
change to this row or its weight needs a new identifying/measurement contract;
it must not simply be dropped because it is difficult to fit.

The new panel, if its launch and smoke succeed, changes one structural
coordinate at a time around the selected point and checks selected responses
at half step. Each point solves prices and the stationary population and
re-normalizes the child-preference scale under the unchanged contract. It can
show which existing parameters move several residuals together, whether local
responses are stable, and where further search would be uninformative without
a specification or measurement change. It cannot validate the diagnostic
income candidate or establish global identification.

## Geography comparison already available

The same saved ACS extraction produces both geographic definitions. These are
alternative measurement systems, not interchangeable targets in a single loss.

| Housing moment | Current 42 metros | National |
|---|---:|---:|
| Mean rooms, capped at nine | 5.561097 | 5.607886 |
| Ownership, heads aged 30–55 | 0.648334 | 0.676260 |
| Recent-parent ownership gap | 0.162896 | 0.127608 |
| Rooms, 3+ versus 1–2 resident children | 0.347067 | 0.385100 |

Source: `../housing_profiles_v1/full/target_recomputed.json`; the current four
rows reproduce exactly. The national comparison has not been substituted into
the calibration or assigned new weights. Its higher ownership target means a
geographic switch is not automatically a remedy for weak ownership fit.

## Work still required before adoption

1. Finish and review the controlled fit panel with every target, actual bounds,
   numerical gates, exact selected repetitions and the unchanged 17 plots.
2. Close the earnings measurement/period/grid/entry contract. If the numerical
   approximation changes, jointly refit the structural parameters under that
   candidate rather than rank specifications at old parameter values.
3. Resolve the target-population description and any necessary remeasurement.
   Keep lifecycle profiles visible as validation; scoring new rows requires
   their identifying role, uncertainty and weights to be specified.
4. Only then decide whether a missing economic margin is needed. A weak local
   fit or a small policy response alone is insufficient reason to add one.

## Numerical execution receipt

Submitted at 23:28 EDT: smoke **18153070**, dependent production **18153071**.
The local three-test suite and independent mocked 28-evaluation loop pass.
The real native smoke and production remain pending review. Replace this
paragraph with terminal outcomes, actual counts and verified fit/sensitivity
findings before the morning report. Submission alone is not numerical success.
