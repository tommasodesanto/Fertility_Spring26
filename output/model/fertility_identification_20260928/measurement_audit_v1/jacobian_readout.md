# What the September 28 Jacobian establishes

Reference: **2007 stationary reference — block0506, September 28 verified
export**. No new calibration is reported here. The common-primary loss remains
19.581. [All 14 fit rows, including weights and contributions, and all 31
parameter estimates/restrictions/bounds](README.md#complete-reference-fit-and-parameter-restrictions)
remain the reference tables.

The Jacobian is the table of changes in measured outcomes when each fitted
parameter changes a little. All 40 actual probes completed: ten coordinates,
positive and negative perturbations, each at two step sizes. Each probe
readjusted the child-benefit level to completed fertility 2.1 and enforced
demographic renewal. These are total calibration responses, including that
adjustment and housing-market clearing, not fixed-benefit causal effects.

## Main fertility finding

Locally, the effective ways to increase the number of children by 25 also tend
to bring first births forward. Mean first-birth age already fits closely:
25.933 in the reference against 25.976 in the data. Thus a move that helps the
early count can damage a target that currently fits.

For example, the full-step derivatives predict that raising the later-birth
choice scale by 10% adds 0.015 children by 25 but lowers first-birth age by
0.199 years. Raising the first-birth scale by 9.865% offsets that timing change,
leaving a gain of only 0.007 children by 25. The half-step derivatives imply a
9.742% offset and a 0.007 gain. The current early-count shortfall is 0.274.
This calculation holds only first-birth mean age fixed to first order; it does
not hold every other target fixed or establish a maximum attainable gain.

The choice scales govern dispersion in the wait/try decisions; they are not
direct child-benefit parameters. First and subsequent births are linked in
household decisions. The separately adjusted child-benefit level also connects
their calibration responses. The Jacobian measures the combined response; it
does not separately assign it to those mechanisms.

## All ten coordinates: illustrative local changes

The entries below are derivative times the stated parameter change, using the
saved central differences. They are linear predictions, not newly solved cases.
Ten-percent changes are larger than the derivative probes. Patience uses 1%
because a 10% increase would exceed its approved bound.

| Parameter increased | Increase | Change in children by 25, full / half step | Change in first-birth age, full / half step |
|---|---:|---:|---:|
| Housing supply level, H0 | 10% | +0.002 / +0.003 | -0.023 / -0.025 years |
| Annual patience factor, beta | 1% | +0.009 / +0.009 | -0.179 / -0.178 years |
| Ownership preference, chi | 10% | -0.005 / -0.005 | +0.056 / +0.057 years |
| First-birth fixed cost | 10% | -0.003 / -0.002 | -0.046 / -0.048 years |
| First-birth choice scale | 10% | -0.008 / -0.008 | +0.202 / +0.202 years |
| Later-birth choice scale | 10% | +0.015 / +0.015 | -0.199 / -0.197 years |
| Bequest strength, theta0 | 10% | +2.598e-06 / +3.339e-06 | -1.453e-04 / -1.503e-04 years |
| First-child housing loading | 10% | -0.001 / -0.001 | +0.009 / +0.009 years |
| Child-benefit curvature | 10% | +1.914e-04 / +1.884e-04 | -4.547e-04 / -4.347e-04 years |
| Tenure-choice scale | 10% | +2.714e-04 / +2.692e-04 | -0.002 / -0.002 years |

The raw [complete derivative table](../lead_review/jacobian.csv) contains all
14 measured outcomes and the derivative of normalized child benefit for every
coordinate. The small response to housing supply does not prove that housing
is economically unimportant: this is one local, normalized equilibrium response.
Similarly, the small relative curvature move at its current low level does not
prove that large absolute curvature changes cannot matter.

## What actual search cases add

The repeated candidate from the experiment fixing the later-birth scale at
twice its reference value, while searching the other nine coordinates, has
early fertility 0.688 but mean first-birth age 23.693. Its common-primary loss
is 807.454. This is a concrete observed trade-off, not a derivative prediction.
It is not the reference, nor proof of the best attainable fit under that fixed
scale. [Complete candidate fit](../resume_v1/selected_export/profile_double/target_fit.csv)
and [all parameters/bounds](../resume_v1/selected_export/profile_double/parameters.csv)
retain the experiment's own weights; the quoted loss rescales every row to
the original primary weights. Primary calibration did not improve in the
completed search. Finite search failures/timeouts do not prove infeasibility.

## Identification and numerical use

The ten-by-ten scored-moment derivative matrix is formally full rank, but
several parameter combinations change the targets very little. The weakest
combinations mostly involve child-benefit curvature, tenure-choice scale and
bequest strength. Their exact directions move with the finite-difference step.
That warns against treating those parameter estimates as precisely identified.
It is not a statistical rank test or proof of global nonidentification.

The full/half weighted condition numbers are 2.695e+06 and 4.654e+05. Seven
half-step singular values exceed the norm of the difference between the two
matrices. The core fertility-direction signs are stable across both steps;
that local trade-off is better supported than precise inversion of the weakest
directions. The older wealth/income dispersion derivative is especially
unstable, but that row has zero calibration weight.

A bounded, damped Gauss–Newton proposal uses this derivative table to choose
a joint parameter move. It predicts reducing the primary loss from 19.581 to
13.335 (13.346 with the half-step derivatives), while early fertility changes
from 0.535 to 0.534. Thus it is a possible improvement to other target fits,
not an early-fertility remedy. It has not been solved or adopted.
[All 14 predicted rows](proposed_step_predictions.csv),
[proposed parameters and bounds](proposed_step_parameters.csv), and the
[bounded unlaunched plan](bounded_experiment_plan.json) are retained.

## What the new spacing diagnostic adds

The Jacobian varied existing parameters while retaining one birth per
four-year period. It did not test that restriction. The separate
[two-birth fixed-policy replay](../two_births_v1/README.md) tests its mechanical
effect at the saved choices. It does not yet solve a new household decision
problem or normalize lifetime fertility back to 2.1. No conclusion about the
2023 transition follows from either exercise.
