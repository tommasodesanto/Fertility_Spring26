# Calibration identification review

Prepared September 29, 2026 for Tommaso De Santo. This reviews existing Jacobian
evidence and recommends changes to the calibration design. An independent Astra
critique was checked against the saved matrices, search source and measurement
receipts. No model was imported or solved; no specification or target was changed.
Calibration launches remain paused during the author-led credit revision.

The principal conclusion is that the current parameters can fit most targeted
averages, but the data do not convincingly distinguish every mechanism. Child-
benefit curvature is the clearest candidate for a parsimonious restriction test.
The age-25 child count should remain informative: replacing it with motherhood
would remove the margin on which the model currently fails.

## Which evidence is being reviewed

The adopted reference remains **2007 stationary reference — block0506,
September 28 verified export**, in
`output/model/fertility_identification_20260928/resume_v1/selected_export/primary/`.
Its common-primary loss is 19.581310760138322. Its checkpoint identity is
`b15ba92dc60e3d5590d2beb6e05d36f71d17b20b1a432edc2c2db926a217309d`.

The later original-model experiment E01 selected `one_birth_024_gn1_0`, with
loss 7.826226594410982. It is not adopted. It retains the original one-birth
architecture and target system, but re-estimates the ten searched coordinates
and re-normalizes the child-benefit level. Its checkpoint identity is
`308245b919c8a50d25972b29e55b5d4fca1ae57d01155ff6050724dc61610451`.
The tables below describe E01, not a newly promoted reference.

| Derivative evidence | Point and method | Finding | Limitation |
|---|---|---|---|
| September 28 reference analysis | R00; 40 probes, central differences at two step sizes | Condition numbers 2.695 million / 465,400; seven half-step singular values exceed the full–half matrix difference norm | That comparison is a numerical sensitivity heuristic, not a statistical rank test |
| E01 overnight final Jacobian | `one_birth_011_gn0_0`; ten one-sided probes, steps 0.05 of declared coordinate scales | Rank 10/10 at relative cutoff (10^{-6}); condition 20,147.54; smallest singular value 0.0046353 | Predates selected E01; derivative noise and step robustness not measured at this point |
| E02 two-birth final Jacobian | `two_birth_012_gn0_1`; same one-sided method | Rank 9/10; smallest singular value 0.00008564 | Smallest value only 12.6% below the numerical cutoff; not proof of statistical nonidentification |

The condition numbers cannot be compared as a measure of improvement: the
points, economic specification in E02, column scaling and derivative method
differ. R00 already has a full/half-step check; the missing check is at the
later selected candidate and, eventually, the revised valid credit specification.

## What the Jacobian means

A Jacobian records how each model moment changes when a parameter changes a
little. Here it also includes housing-market adjustment and the separate outer
loop that adjusts the child-benefit level to completed fertility 2.1.

Let $\theta$ collect the ten searched parameters, let $\psi$ be the normalized
child-benefit level, and let $F(\theta,\psi)=2.1$ denote the normalization.
Provided $F_\psi\ne 0$ and local differentiability holds, the measured derivative is

$$
J=W^{1/2}\left[m_\theta-m_\psi\frac{F_\theta}{F_\psi}\right]S.
$$

Here $m$ collects scored moments, $W$ contains their working weights, and
$S$ contains the declared parameter scales. The derivatives of $m$ and $F$
already include equilibrium adjustment. A column therefore describes a move
along the fertility normalization, not a causal response holding $\psi$ fixed.
This explains why an increase in a birth-choice scale can affect first births,
later births and the normalized benefit simultaneously. A fall in normalized
$\psi$ alone does not establish a lower structural valuation of children.

Ten scored moments for ten searched coordinates, plus a separate normalization
for $\psi$, satisfy the counting requirement. They do not establish that the
moments contain enough independent information, nor that the working weights
are an efficient SMM covariance matrix. We have no statistical confidence
intervals for these fitted parameters from this analysis.

In E01, the weakest scaled direction has coefficients 0.9893 on curvature,
0.1039 on tenure-choice dispersion and 0.0963 on the first-birth cost. The signs
of a singular vector are arbitrary. These loadings describe a combination that
barely changes weighted moments near the round center; they are not parameter
standard errors. Column scales are hand-selected, including 0.1 for curvature
and 0.01 for annual patience. Report robustness to scales rather than calling
curvature globally unidentified.

## Are the parameters the right ones

The following mappings are economic rationales, not one-moment identification
proofs. All parameters affect several moments through household choices and
equilibrium.

| Parameter | Intended information | Judgment before further estimation |
|---|---|---|
| Housing supply level $H_0$ | Mean rooms and the equilibrium housing quantity | Retain; an average quantity is useful, but the inferred level is conditional on housing-demand primitives and the fixed supply elasticity |
| Annual patience $\beta$ | Wealth relative to earnings | Retain; inspect wealth by age because an aggregate fit can conceal incorrect saving paths and credit constraints |
| Owner housing-service premium $\chi$ | Ownership at ages 30–55 | Retain provisionally; credit access and tenure dispersion can substitute for an ownership preference |
| First-birth fixed utility cost | First-birth timing, with childlessness and young-age fertility | A distinct entry-to-parenthood preference, but aggregate moments may confound it with taste dispersion; the zero-cost test does not prove it indispensable |
| First-birth taste scale | Childlessness and the distribution of first-birth timing | Retain provisionally; distinguish dispersion from a permanent preference for motherhood and from omitted conception risk |
| Later-birth taste scale | One-child share among mothers and additional-birth timing | Retain provisionally; birth-order-specific timing and spacing are more informative than one terminal count share alone |
| Bequest strength $\theta_0$ | Annual child-directed bequest flow relative to wealth | Retain, but qualify the mapping: the flow also depends on wealth, mortality and who has children; check late-life wealth and ownership |
| First-child housing-share loading | Rooms around first birth | Retain provisionally; the active target is the 1.465-room matched ACS proxy, not the older 0.720 PSID event-study estimate; validate the observer and empirical comparison |
| Child-benefit curvature | Young child counts and the distribution across numbers of children | First candidate for a restriction/profile test; it dominates the weakest E01 direction and its economic object is concurrent children at home |
| Tenure-choice taste scale | Ownership differences around parenthood, jointly with the ownership level | A response/dispersion parameter needs gradients or response variation; an aggregate ownership level alone does not establish its value |

The child-benefit specification is $B(m)=\psi m^{1-\gamma}$, where $m$
is children currently at home. Curvature governs the benefit of concurrent
children, while the empirical fertility stock counts children ever born.
Birth spacing and children leaving home connect these two objects. Thus
curvature should be assessed together with the child-departure process,
not treated as a direct parameter of the observed completed-child distribution.
Its marginal benefits are $B(2)-B(1)$ and $B(3)-B(2)$; report these alongside
the normalized benefit level in any future profile.

For example, E01's normalized level 0.121845517 and curvature 0.065734962
imply direct per-period utility increments of 0.12185, 0.11099 and 0.10723
for the first, second and third child concurrently at home. The benefit is
already close to linear. These are utility increments, not cash amounts or
net incentives after material costs and continuation values. Changing curvature
also changes the normalized level, so this arithmetic is not a refitted result.

The age-25 count has a legitimate literature precedent. Sommer (2016, section
4.5) targets a mean birth count of 0.8 at age 25 to discipline a curvature
parameter governing the fertility age profile. Her model includes child-quality
production and time/money inputs. That supports the moment's rationale, but
does not transfer identification to our curvature in benefits from children
at home. [Author-hosted paper](https://www.kamilasommer.net/Fertility.pdf).

The first-birth cost shifts the value of becoming a parent. Taste scales govern
dispersion in wait/attempt decisions and also affect continuation values through
the log-sum term. They are therefore more than harmless numerical smoothing.
Removing them, equating first/later scales, or introducing permanent fertility
tastes would each be an economic specification change. None is adopted here.

The exact-zero-cost E03 experiment raised loss from 7.826 to 518.306 when other
coordinates were held and $\psi$ re-normalized. A bounded nine-parameter local
refit reached 37.076. This establishes deterioration in the tested deletion
and local refit, not global necessity of the cost. Neither a normalization
censor nor an unsuccessful bounded refit proves infeasibility.

## What the fertility miss tells us

| Matched age-25 object | Data | E01 |
|---|---:|---:|
| Children ever born per woman, capped at three | 0.810 | 0.530 |
| Share who are mothers | 45.7% | 44.8% |
| Children per mother, capped at three | 1.770 | 1.185 |

The extensive margin is close; the conditional count is the main miss. At
matched age 26 the count is 0.923 in data and 0.609 in E01, so an age relabeling
does not resolve it. Ages 18–19 are modeled; pre-18 births are omitted.
CPS age-17 women average 0.084 children, and NCHS period first births below
18 comprise 7.731% of all first births. These are different cohort stocks and
period flows. Their contribution to the age-25 gap, including subsequent births
to early mothers, has not been established.

E02 permits an additional birth opportunity within a model cell, but adds an
independent taste opportunity and uses a common event-time proxy. Its selected
age-25 count is 0.606; conditional children improve to 1.385 while motherhood
falls to 43.8%. Overall loss is 7.842. This weakens the claim that one birth per
cell is the sole explanation; it does not certify a clean spacing mechanism
or integration into a calendar-time transition.

Claude's 0.665 ceiling conditions on empirical first-birth cell shares and
additional cohort assumptions. Matching mean first-birth age does not fix those
shares. Neither that calculation nor the local Jacobian proves the active
target system globally infeasible. The mean leaves considerable freedom in the
shape of the timing distribution.

The stationary baseline deliberately approximates 2007 with replacement
fertility. The period timing rows, older CPS child-count rows and younger CPS
age-25 stock are not one cohort's history. Quantify this approximation instead
of silently substituting a different fertility level or target. Nor should
terminal normalization be equated arithmetically with the capped age-40–44
stock. A later 2023 transition fit is a separate test; it cannot certify this
baseline Jacobian in advance.

## Recommended sequence after the code review

Working weights also need interpretation. The age-25 row uses weight 100,
equivalent to a working scale of 0.1 children, whereas its person-bootstrap
standard error is about 0.028. This does not establish that inverse-variance
weighting is appropriate: the bootstrap is not a survey-design variance, and
stationary approximation and observer error are separate concerns. Document
these choices and estimate joint empirical uncertainty before presenting
precision or an overidentification test. Changing weights can move the fitted
trade-off; it does not create identifying information.

1. **Establish one valid economic specification.** Finish the author-led credit
   revision, reconcile the retained negative-debt entrants, and authenticate
   code, entry distributions, saved parameters and observers. Preserve the old
   reference. Existing Jacobians remain evidence about the old specification.

2. **Test parsimony before adding parameters.** Profile several explicit fixed
   curvature values, including zero, re-estimating the remaining coordinates
   and retaining all ten scored moments, three validation rows, fertility 2.1
   normalization and renewal. Nine free parameters with ten scored moments
   avoid a counting deficit. Judge complete fit and economically important
   predictions, not loss alone. This is a proposal, not authorization to run.
   The separate first/later zero-taste tests remain paused; they should not be
   bundled with credit or curvature changes without a disclosed design.

3. **Improve the information before replacing targets.** First measure the
   distribution of first-birth ages, young child-count distributions and
   birth-order-specific spacing on compatible cohorts. Show pre-entry births
   separately. Assess fertility by pre-birth resources or education before
   adding a mechanism to explain a suspected income gradient. Keep age-25
   motherhood and conditional children as diagnostics initially. Their product
   equals children per woman, so adding all three does not create three
   independent identifying facts. A revised objective requires their joint
   uncertainty and an explicit target contract.

4. **Certify local derivatives at the selected valid point.** Repeat the center,
   use central differences at two step sizes where feasible, report boundary
   one-sided derivatives separately, and compare step/repeat variation with
   the weakest sensitivities. A complete two-step central check requires 40
   probes plus center repeats; it is a separate bounded Torch diagnostic,
   not a search restart. Solve small displacements along the weak direction
   before trusting its local inversion. All future budgets need fresh timing
   evidence and explicit limits.

5. **Use derivatives to improve search, not to certify economics.** A bounded
   trust-region or damped Gauss–Newton method is preferable to blindly widening
   random search when derivatives are reliable. Record predicted and realized
   changes in every moment, and shrink steps after poor prediction. An SVD
   ridge limits moves in weak directions; it does not identify their parameters.
   Derivative-based $\psi$ warm starts must still reach the same normalization
   root and pass the same gates. Speed gains remain a benchmark question.

My priority is steps 1–3. The point of better numerics is to make the chosen
model inspectable and estimable; it cannot replace better identifying evidence.

## Complete E01 fit and searched parameters

Values below are rounded for readability. The exact complete 14-row comparison
is in [full target fits](../../output/model/fertility_identification_20260928/zero_first_birth_cost_v1/readout_v1/full_target_fit.csv),
using its **Original selected** columns. The original scored loss is
7.826226594410982, of which 7.788636564687688 comes from the age-25 count.

| Scored moment | Target | Model | Model minus target | Weight | Loss contribution |
|---|---:|---:|---:|---:|---:|
| Childlessness, women 40–44 | 0.198279 | 0.198177 | -0.000102 | 35532.304 | 0.000367 |
| One child among mothers 40–44 | 0.213655 | 0.214028 | +0.000373 | 26952.821 | 0.003753 |
| Mean first-birth age | 25.976264 | 25.965339 | -0.010925 | 139.828 | 0.016689 |
| Wealth / earnings | 6.926584 | 6.896046 | -0.030537 | 7.595 | 0.007083 |
| Annual child-directed bequests / wealth | 0.00729102 | 0.00728053 | -0.00001049 | 5165289.256 | 0.000568 |
| Mean rooms | 5.729434 | 5.729955 | +0.000521 | 128.021 | 0.000035 |
| Ownership, ages 30–55 | 0.676260 | 0.676202 | -0.000058 | 2339.362 | 0.000008 |
| First-birth rooms difference | 1.465000 | 1.463720 | -0.001280 | 137.565 | 0.000225 |
| Recent-parent ownership difference | 0.127608 | 0.127036 | -0.000572 | 27055.823 | 0.008862 |
| Children per woman, age 25 | 0.809528 | 0.530446 | -0.279081 | 100.000 | 7.788637 |

| Untargeted validation moment | Data | Model | Model minus data | Weight | Loss contribution |
|---|---:|---:|---:|---:|---:|
| First births at age 30 or later | 0.249278 | 0.225697 | -0.023581 | 0 | 0 |
| Older wealth/income dispersion | 3.515935 | 2.919479 | -0.596456 | 0 | 0 |
| Family rooms difference | 0.385100 | 0.473121 | +0.088021 | 0 | 0 |

Completed fertility **(normalization)** is 2.099609 against 2.1, gap -0.000391;
it has no scored weight or loss contribution and passes its tolerance.

| Searched parameter | E01 estimate | Original bounds | Bound assessment |
|---|---:|---|---|
| $H_0$ | 6.104357 | [0.2, 80] | Interior |
| Annual $\beta$ | 0.969126 | [0.94, 0.99] | Interior |
| $\chi$ | 1.098496 | [0.1, 5] | Interior |
| First-birth fixed cost | 0.352709 | [0, 8] | Interior |
| First-birth taste scale | 0.108549 | [0.02, 50] | Not binding; generic wide-interval near-bound flag |
| Later-birth taste scale | 0.222398 | [0.02, 50] | Not binding; generic wide-interval near-bound flag |
| $\theta_0$ | 0.156232 | [0, 8] | Interior |
| First-child housing loading | 0.115404 | [0, 0.25] | Interior |
| Child-benefit curvature | 0.065735 | [0, 0.8] | Interior |
| Tenure-choice taste scale | 0.012684 | [0.001, 0.1] | Interior |

The normalized $\psi$ is 0.121845517; positivity is required and it is not an
eleventh unrestricted search coordinate. The taste-scale flags use 1% of a
very wide interval and should not be interpreted as estimates approximately
zero. [All 31 parameter estimates, fixed restrictions and original flags](../../output/model/fertility_identification_20260928/zero_first_birth_cost_v1/readout_v1/full_parameters.csv)
remain available. Both final E01 repeats passed. The 17 standard diagnostic
plots remain retained; this review neither replaces nor changes them.

## Evidence locations

- [Reference central-difference readout](../../output/model/fertility_identification_20260928/measurement_audit_v1/jacobian_readout.md)
- [Raw reference derivatives](../../output/model/fertility_identification_20260928/lead_review/jacobian.csv)
- [Original overnight round-center matrix](../../output/model/fertility_identification_20260928/two_stream_overnight_v1/morning_readout_v1/run_v1/one_birth/jacobian_1.json)
- [Declared scaling and search algorithm](../../output/model/fertility_identification_20260928/two_stream_overnight_v1/search.py)
- [Round-center weak directions](../../output/model/fertility_identification_20260928/two_stream_overnight_v1/comparison_v1/weak_direction_round_centers.csv)
- [Age and teenage-birth measurement](../../output/model/fertility_identification_20260928/measurement_audit_v1/age_tail_v1/README.md)
- [Child-benefit implementation](../../code/model/intergen_eqscale_seq_optimized/child_preferences.py)
- [Experiment register](../../output/model/fertility_identification_20260928/EXPERIMENT_REGISTER.md)

The current credit decision and remaining implementation/entry gaps belong to
`CALIBRATION_STATUS.md`. That live status takes precedence over historical
calibration contracts. No new numerical evidence about the revised credit model
was generated for this review.
