# Audit against the circulated quantitative model

Date: 2026-08-05

## Bottom line

Yes. We have been cutting corners in the way results were promoted and
reported.

The most important error was to call E5b the maintained model without saying
that it is **not the model in the July circulation**. The circulated paper uses
Stone--Geary needs and a one-shot completed-family-size choice. E5b is a later,
explicitly provisional fork with an imposed equivalence scale and sequential
birth attempts. The circulated code was not overwritten, but E5b should have
been called an experimental specification, not the paper benchmark.

The E5b policy packet is still further from the circulated contract. It fixes
adult population mass, injects entrants exogenously, does not close the
government budget, and does not compute welfare. It is a mechanism diagnostic,
not the population-adjusted policy experiment described in the paper.

Returning mechanically to the old circulated estimates is also unsafe. The
circulated M5 wealth readout predates the July 23 balance-sheet-timing repair,
its entry elasticity was imposed rather than estimated, and the draft contains
two different sets of policy numbers. The correct next benchmark is therefore
a **circulation-consistent, repaired and recertified M5**, not either the saved
old table or E5b.

## Audit standard

The comparison object is the July circulated paper and slides:

- `latex/intergenerational_housing_fertility_paper_draft.tex`
- `latex/intergenerational_housing_fertility_note_slides.tex`
- `latex/quantification_rewrite_review.tex`
- the entry and population equations in `latex/model_writeup.tex`

The audited quantitative object is the certified E5b estimate and its policy
packet:

- `output/model/eqscale_seq_e5b_recalibration_20260725/`
- `output/model/eqscale_seq_e5b_policy_packet_20260725/`
- `code/model/intergen_eqscale_seq_optimized/`

Each item below is classified as an exact match, a disclosed approximation, an
undisclosed or insufficiently disclosed shortcut, or a contradiction.

## Priority findings

| Priority | Object | Finding | Classification | Consequence |
|---|---|---|---|---|
| 1 | Model identity | Circulated M5 is Stone--Geary with a one-shot family-size choice; E5b uses equivalence-scale preferences and sequential birth attempts. | **Contradiction if E5b is called the paper benchmark.** The changes were approved as experiments and documented in the July 24 provisional note. | E5b evidence cannot be inserted into the circulated paper as a recalibration of the same model. |
| 1 | Population and entry | The paper has outside-option entry and endogenous stationary city scale. E5b injects `1/J` entrants and renormalizes population to one after every forward pass. Mature local children are measured but do not determine the next entrant mass. | **Contradiction.** | E5b policy results contain composition and price responses only; they are not long-run population-adjusted effects or changes in total city births. |
| 1 | Fiscal closure | The circulated equilibrium requires government accounting. E5b pays the purchase grant without a financing source and discards property-tax revenue. | **Contradiction.** | The grant and tax-plus-grant cases are not feasible government policies and cannot support welfare or incidence claims. |
| 1 | Circulated M5 measurement | The saved M5 wealth moments use the pre-repair timing convention that could pair inherited liquid wealth with newly chosen tenure before the housing transaction. | **Known unresolved defect in the circulated benchmark.** | The circulated fit cannot simply be restored; its wealth block must be remeasured and the model recalibrated. |
| 1 | Reported policy numbers | The circulated draft reports `+2.1%/+2.7%` total-birth effects in its table, but later in the same section and conclusion reports `+3.4%/+4.9%`. | **Internal contradiction in the circulated document.** | There is no unique circulated policy result to restore without tracing the exact source records. |
| 2 | Sequential child maturation | The absence of separate child ages is intentional and consistent with the current specification. The defect is narrower: E5b puts all children on one shared maturation clock, so they leave home as a group. | **Implementation shortcut inconsistent with the author’s intended specification.** | The model must retain parity as children ever born but stochastically mature each child independently from the existing count-at-home state. |
| 2 | Grant eligibility | The E5b grant is available whenever a renter with children at home buys an eligible home. This matches the circulated eligibility as written, but there is no state recording prior receipt. | **Exact as written; unresolved design ambiguity.** | A household can in principle receive the grant again after returning to renting. If the intended policy is first-time or once-per-family, the current code cannot implement it. |
| 2 | Identification | E5b has 12 moments and 10 free parameters, but the gate checks only the count. No fresh full 12-by-10 E5b Jacobian or weak-direction audit is in the certified packet. | **Incomplete validation.** | Formal overidentification by counting does not establish that the ten parameters are separately informative. |
| 2 | Objective weights | Six E5b rows use measured or declared standard errors; the others use a synthetic standard error equal to 5% of the target. The two timing rows use declared window-spread values, not a joint sampling covariance matrix. | **Disclosed approximation, previously overinterpreted.** | Loss `385.14` is a calibration criterion, not a formal SMM $J$-statistic. |
| 2 | Tenure choice | The circulated M5 estimates a positive tenure-choice scale near `0.01`. E5b fixes it to zero, making tenure deterministic conditional on the state. | **Author-approved experimental restriction, but a contradiction with the circulated specification.** | Purchase responses occur at sharp wealth thresholds and are especially sensitive to grids and the discrete house menu. |
| 2 | Housing grid | E5b retains five owner sizes, `[2,4,6,8,10]`. No E5b dense-menu recalibration or policy robustness exists. The earlier M5 fixed-parameter audit changed ownership from `0.658` to `0.713` when the menu was densified. | **Incomplete numerical robustness.** | The very large ownership response to the E5b grant may partly reflect discrete threshold geometry. |
| 2 | Wealth grid | E5b is estimated at `Nb=120`; I found no `Nb=240` winner verification, despite the standing search-at-120/verify-at-240 convention and earlier material grid drift. | **Incomplete numerical robustness.** | Threshold, ownership, and wealth-tail moments are not yet certified as grid-stable. |
| 2 | Entry elasticity in circulated policy | The later M5 entry/scale implementation uses benchmark entry probability `0.5`, entry taste scale `2`, and local-birth weight `1`. These objects are not disciplined by the calibration. | **Disclosed but economically undisciplined closure.** | Even the population-adjusted M5 policy results are sensitivity exercises; the roughly 29% scale response in the funded rerun is not paper-ready. |
| 3 | Timing measurement | The model has four-year age points but is fitted to continuous NCHS mean-age and age-30-share moments. | **Approximation.** | The timing fit is partly a binning comparison; the declared timing uncertainty does not replace a full interval-mapping exercise. |
| 3 | Family ownership gap | The empirical row compares households by observed children, while the model comparison is new parents versus nonparents under its state definitions. | **Disclosed near-match, not identity.** | This row is useful but should not be described as exact measurement equivalence. |
| 3 | Cohort alignment | Fertility uses 1979--84 birth cohorts, housing uses mainly 2019--23 cross-sections, and wealth uses a longer PSID window. | **Disclosed maintained assumption.** | The calibration assumes stable preferences and interpretable cross-cohort moments; this requires a robustness statement, not concealment. |

## What is implemented faithfully

The audit did not find evidence that all model work is corrupted. Several
central objects are correctly implemented:

1. The E5b equivalence-scale equation matches its provisional mathematical
   note. With $e(n)=((2+0.7n)/2)^{0.7}$, CRRA utility uses the correct
   multiplier $e(n)^{\sigma-1}$, and the expenditure share is
   \[
   \alpha(n)=\alpha_0-\delta_{\mathrm{jump}}\mathbf 1\{n\ge1\}-\delta_a n.
   \]
2. The literal parity states $0,1,2,3+$, the no-two-births-in-one-period
   rule, fecundity splitting, and mass conservation have direct tests.
3. The home-purchase feasibility condition uses the correct cash requirement
   \[
   b\ge (1-\phi)pH,
   \]
   with $\phi$ as the financed share.
4. The entrant wealth distribution, post-retirement survival schedule, sale
   cost, rental cap, owner menu, and upward-sloping housing supply are active.
5. The E5b winner reproduces all twelve stored model moments and its scalar
   loss exactly under repeated strict solves.
6. The rejected E6 fecundity, permanent-income, and readiness mechanisms are
   default-off. They did not overwrite the pre-E6 behavior.
7. The circulated M5 code remains in the repository. The problem is promotion
   and validation, not loss of the old model.

## E5b calibration evidence

E5b's certified loss is `385.142880159`, with market residual
`1.109e-6`. Eight rows are close, but childlessness, late first births,
aggregate wealth, and late-life wealth dispersion fail materially.

| Moment | Target | Model | Gap | Weight | Loss contribution |
|---|---:|---:|---:|---:|---:|
| Completed fertility | 1.918000 | 1.903573 | -0.014427 | 1,425.739 | 0.297 |
| Childlessness | 0.188000 | 0.080188 | -0.107812 | 17,180.744 | 199.700 |
| Mean age at first birth | 25.310561 | 25.637433 | 0.326872 | 16.000 | 1.710 |
| First births at age 30+ | 0.270062 | 0.332334 | 0.062272 | 15,625.000 | 60.590 |
| Rooms added at first child | 0.664435 | 0.657205 | -0.007230 | 906.056 | 0.047 |
| Rooms gap, 3+ versus 1--2 children | 0.367700 | 0.372240 | 0.004540 | 2,958.515 | 0.061 |
| Family ownership gap | 0.167662 | 0.165878 | -0.001783 | 14,229.591 | 0.045 |
| Ownership rate, ages 30--55 | 0.575472 | 0.563758 | -0.011714 | 1,207.846 | 0.166 |
| Mean occupied rooms | 5.779970 | 6.075505 | 0.295534 | 11.973 | 1.046 |
| Wealth / annual gross labor earnings | 6.873100 | 5.454754 | -1.418346 | 6.288 | 12.649 |
| Annual bequest flow / wealth | 0.008800 | 0.009078 | 0.000278 | 5,165,289.256 | 0.399 |
| Old wealth p90/p50 | 3.448111 | 2.068370 | -1.379741 | 56.960 | 108.433 |

Full machine-readable table:
`output/model/eqscale_seq_e5b_recalibration_20260725/report/target_fit_full.csv`.

### Free parameters and restrictions

| Parameter | Estimate | Search bounds | Bound flag |
|---|---:|---:|---|
| Annual discount factor | 0.991097 | [0.9400, 0.9995] | No |
| Per-child expenditure-share shift | 0.012163 | [0, 0.25] | No |
| First-child expenditure-share jump | 0.023896 | [0, 0.25] | No |
| Child utility flow | 0.585018 | [-3, 3] | No |
| First-birth taste scale | 3.552006 | [0.02, 50] | No |
| Continuation taste scale | 1.254006 | [0.02, 50] | No |
| Ownership utility premium | 1.038033 | [0.10, 5] | No |
| Housing-supply scale | 8.162391 | [0.20, 80] | No |
| Bequest strength | 0.168435 | [0, 8] | No |
| Bequest curvature shift | 0.050046 | [0.02, 16] | **Near lower bound** |

Fixed outside the E5b estimation are the consumption share
`alpha_cons=0.733`, deterministic tenure `tenure_choice_kappa=0`, and
child-invariant bequest utility `theta_n=0`.

Full machine-readable table:
`output/model/eqscale_seq_e5b_recalibration_20260725/report/parameter_bounds.csv`.

## Policy audit in equations

The circulated counterfactual requires

\[
Q(S,p)=M+S B_0(p), \qquad E_i=Q(S,p)\pi_i^E(p),
\]

and stationary scale

\[
S(p)=\frac{q^E(p)M}{E_0(p)-q^E(p)B_0(p)}.
\]

The benchmark sets $S^*=1$, calibrates the outside value to a benchmark
entry rate, recovers $M$, and then holds $M$, the outside value, and the
entry-shock scale fixed in counterfactuals.

E5b instead uses

\[
E=1/J, \qquad G\leftarrow G/\int dG,
\]

after the forward pass. That is a fixed-population normalization. The code
computes mature-child flows but does not feed them into $E$.

The fiscal policy should also satisfy, for an equal rebate $T$,

\[
\tau^p p\int H^O(x)\,dG(x)
=T\int dG(x)+A\Pr(\text{eligible purchase})\int dG(x).
\]

The E5b packet solves neither $T$ nor this identity. Consequently its
reported policy effects are correctly interpreted only as unfunded,
fixed-population comparative statics:

| Case | Completed fertility | Change | Ownership | Change | Price |
|---|---:|---:|---:|---:|---:|
| Baseline | 1.903573 | -- | 0.563758 | -- | 0.815625 |
| Grant | 1.905159 | +0.001586 | 0.635743 | +0.071985 | 0.819360 |
| Tax plus grant | 1.906023 | +0.002450 | 0.619109 | +0.055351 | 0.661504 |

The response is almost entirely tenure, not fertility. There is also no
tax-only E5b arm, so the package cannot separate the tax effect from the grant
effect.

## Required correction before further policy work

1. Treat the July circulated M5 as the paper's economic specification and
   E5b as an experimental alternative.
2. Repair the experimental sequential specification without adding child-age
   states: births map $(n,m)$ to $(n+1,m+1)$ and each of the $m$ children at
   home matures independently with four-year probability $2/9$.
3. Re-estimate under coherent targets, then run a fresh local Jacobian,
   `Nb=240` verification, and dense-owner-menu robustness.
4. Freeze and identify the outside-entry objects or report a predeclared
   sensitivity grid. Do not present an assumed entry elasticity as estimated.
5. Solve a fiscally closed baseline and the tax-only, grant-only, and combined
   policies under the full entry/scale equilibrium. Fixed-population E5b rows
   remain decompositions only.
6. Add a prior-receipt state if the grant is intended to be one-time or for
   first-time buyers.
7. Reconcile the contradictory policy numbers and regenerate every table and
   figure from one certified source record before the paper or slides are used
   again.

Until these steps are complete, no E5b policy elasticity and no saved
population-adjusted M5 total-birth number should be described as a paper
result.

## August 6 follow-up: repaired experimental E5

The two author-directed implementation defects in the experimental E5 branch
have now been repaired and tested. Children mature independently from the
existing count-at-home state, and the policy driver now imposes the government
budget and the circulated outside-option entry/stationary-scale equations. No
child-age state or new estimated parameter was added.

The repaired ten-parameter, twelve-target recalibration is numerically
certified: all eight chains have exact strict repeats, the selected loss is
`249.186326675`, and an independent solve reproduces all moments. This is an
improvement over historical E5b loss `385.142880159`, but childlessness,
completed fertility, and old-age wealth dispersion remain materially wrong.

The funded policy audit also passes. Relative to the rebated 1% baseline, the
2% tax-only case changes completed fertility by `+0.068%`, house price by
`-11.84%`, stationary scale by `+21.06%`, and total births by `+21.13%`.
Adding the `0.4` purchase grant changes those objects by `+0.383%`, `-12.22%`,
`+18.71%`, and `+19.13%`, respectively. These scale results inherit the
externally imposed entry probability `0.5` and taste scale `2`; they are
sensitivity results, not estimated elasticities.

This closes the two implementation items but does not make E5 the circulated
model. The remaining audit requirements are a circulation-consistent repaired
M5 benchmark, identified entry behavior, welfare, a decision on one-time grant
receipt, fresh Jacobian and grid/menu robustness, and reconciliation of all
paper tables. The repaired E5 result is therefore retained as a certified
experiment and is not promoted.
