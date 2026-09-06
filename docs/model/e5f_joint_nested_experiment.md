# Experimental simultaneous nested choice

Author authorization: September 6, 2026, 20:17 UTC. This experiment is the
first quantitative priority. Production is unchanged. The author approves
experimentation, not adoption of a new model or nesting restriction.

## First-stage economic contract

This is a **one-date, fixed-price, fixed-continuation choice diagnostic** on
the reproduced retained September 4 calibration's 2023 pre-choice population.
It is not an equilibrium, calibration, historical housing-response estimate,
policy forecast, or welfare comparison. No twelve-moment loss is computed for
the new choice rule. Every original target and parameter remains recorded with
the benchmark; the new shock scales are diagnostic values, not estimates.

For state x, define q(tau,x) as the maximum deterministic conditional value
over renter or owner housing products in tenure tau, with the original
transactions, budgets, borrowing limits, conditional consumption/saving and
baseline continuation values. Renter rooms are continuous, owner products are
the original five sizes. Product-specific tenure shocks are removed in this
experimental block. The four alternatives are (tenure, wait/attempt).
The retained benchmark has exactly one location. There is no spatial choice
in this experiment; q includes the original one-location amenity/staying
constant. The driver rejects multi-location or readiness specifications.

Housing tenure is selected with the birth attempt. If the attempt succeeds,
housing size, consumption and saving may adjust within that tenure. The same
tenure is maintained on failure and success. This excludes success-contingent
switches between rent and own and must be disclosed separately from revelation
of tastes. It is a coherent first experiment, not the only simultaneous model.
An owner may select different original owner products on success and failure.

Let pi be the original age-specific conception probability; x+ adds one birth
and child at home, and F is the original first-birth utility cost, zero for
continuation births. The four deterministic plan values are

\[
Q_{\tau,0}(x)=q(\tau,x),\qquad
Q_{\tau,1}(x)=(1-\pi)q(\tau,x)+\pi[q(\tau,x^+)-F].
\]

Plans with a positive-probability infeasible outcome are infeasible. Zero-weight
outcomes do not remove a plan. Attempt is unavailable outside fertile ages and
at the top parity; no second birth is applied to mass just born in this date.

All four taste components are jointly revealed before the plan is chosen.
Use the standard nested-GEV law with tenure nests, marginal scale kappa>0 and
within-tenure dissimilarity 0<lambda<=1, so kappa_inner=lambda*kappa:

\[
F(\epsilon)=\exp\left\{-\sum_\tau
  \left[\sum_a e^{-\epsilon_{\tau a}/(\lambda\kappa)}\right]^\lambda\right\}.
\]

The displayed law uses uncentered marginals. Subtract the same Euler-constant
times kappa from all components for mean-zero marginals; choices are unchanged.
The corresponding mean-zero-shock expected maximum and probability are

\[
S_\tau=\lambda\kappa\log\sum_a e^{Q_{\tau a}/(\lambda\kappa)},\quad
V=\kappa\log\sum_\tau e^{S_\tau/\kappa},\quad
p_{\tau a}=\frac{e^{S_\tau/\kappa}}{\sum_s e^{S_s/\kappa}}
\frac{e^{Q_{\tau a}/(\lambda\kappa)}}{\sum_b e^{Q_{\tau b}/(\lambda\kappa)}}.
\]

Nesting is correlation, not an order of economic decisions. Within-tenure
probabilities become deterministic as lambda approaches zero, giving a
conceptual bridge to the simplified theory's conditional renter/owner problems.
This does not equate their distinct fertility and lifecycle architectures.
Reference: Train (2009), chapter 4, pp. 79-83:
https://eml.berkeley.edu/books/choice2nd/Ch04_p76-96.pdf

## Controls, outputs and interpretation

1. Replay the original fixed-price Bellman/forward calculation and require
   its values, conditional policies and distributions to reproduce the saved
   reference. Collect q with observational wrappers only; scientific files stay
   hash-matched. Independently check tenure-specific maxima reconstruct the
   original deterministic housing maximum. Replay with frozen continuation V
   as well: the benchmark's backward recursion must reproduce itself.
2. Apply the joint operator on the exact pre-choice population. Scatter each
   origin's joint plan and conception outcome directly into the existing
   transaction maps; never multiply unrelated unconditional fertility and
   tenure marginals or let incoming birth mass attempt again.
3. Grid: kappa in {.005,.05,.5,2.5}, lambda in {.05,.5,1}; twelve cases.
   Add a lambda=1e-6 case at kappa=.5 for the conditional deterministic limit.
   All scales apply to every age and parity as a diagnostic restriction.
4. At each grid point also evaluate a fertility-first sequential-logit
   construction on the **same four Q values**, using tenure scale kappa and
   fertility scale lambda*kappa. This is a mathematical control with the same
   payoff/commitment menu. It changes the random-utility law along with
   revelation; do not call it a pure timing causal decomposition.
   Its formulas are
   \(\widetilde S_a=\kappa\log\sum_\tau\exp(Q_{\tau a}/\kappa)\) and
   \(\widetilde V=\lambda\kappa\log\sum_a\exp(\widetilde S_a/(\lambda\kappa))\).
   These are sequential expected values, not a fertility-nest GEV law when
   \(\lambda<1\). Probabilities are the outer action probability times the
   conditional tenure probability.
5. Report current birth flows, ownership, housing demand, the fixed-price
   demand/supply gap, attempt/tenure cross-tabulations, by-age and boundary
   exposure, expected-value screens and full readable checkpoint receipts.
   Preserve the reference's seventeen standard graphs and add explicitly
   supplemental joint-choice plots. The new state is not a full model solution.

## Verification and budgets

Before launch: independent arithmetic/finite-difference checks of logsum
probabilities; lambda=1 flat logit; near-zero lambda deterministic conditional
choice; translations; infeasible alternatives; pi=0/1 and mixed conception;
direct four-alternative enumeration and population mass/birth conservation.

One exact-loop smoke first: two independently replayed captures and two grid
cases through the same collection/plot/checkpoint loop. Original grid is 120
wealth nodes, 17 ages, 15 income states, six housing products, four parity
states and four at-home-child states. Each capture makes one frozen-continuation
Bellman call and requires exact equality with the saved ordinary replay.
Budget two calls for the two smoke captures, conservatively 2-5 minutes each including
first JIT; allow 35 minutes total on Torch, one CPU/16GB. Stop immediately on
source, replay, probability, exposure/feasibility or accounting failure.

After the smoke, reuse its saved choice-value capture for the thirteen-case
panel; it requires no additional Bellman call and at most 20 minutes on one
CPU/16GB. No calibration or market search. Write a heartbeat at least every
60 seconds and latest-completed/best-so-far files after each case, with
best-so-far explicitly labeled no calibration selection. Absolute first-stage
cluster compute cap 55 minutes plus queueing. No automatic wider sweep.
An independent reviewer gets 20 minutes; lead implementation/review budget
90 minutes before a concrete progress decision. Failures remain visible.

### Exact-loop smoke finding and bounded revision

The first smoke (Torch 17067221, 43 seconds) exactly reproduced all twelve
reference arrays twice, with identical conditional-value and product hashes.
It stopped on an absolute budget-excess-mass gate at 3.29036219045e-9. All
seven occupied exceptions, including each state's mass and expenditure gap,
are exactly identical to the September 5 benchmark audit. They are the
already documented reporting-floor issue, not new joint-choice exposure.

The revised experiment records the full absolute budget audit and also audits
the positive part of the state-by-state change in destination mass. Its gate
is additional budget-violating mass <=2e-10, preventing offsetting increases
and decreases from concealing a new problem. No household rule or production
gate changes. A new hash contract and snapshot rerun the complete smoke;
the failed first run is retained. Total first-stage cluster time remains
bounded by 55 minutes, including the first 43 seconds.

## Completed first-stage result

Final smoke 17068184 and panel 17068310 completed with exit zero. The final
contract additionally corrects reference birth reporting to distinguish
explicit births from the topcode-adjusted child units also reported for every
experimental case. It changes no choice equation; the exact-loop smoke was
rerun after the correction. The final contract is `final_contract.json` in
`output/model/e5f_joint_nested_experiment_20260906a/`.

All 26 cases pass the stated checks; all four smoke cases reproduce exactly
in the full panel. The maximum joint/control difference is 0.000106336 birth
units per 100 households over four years and 0.0000053473 ownership percentage
points. Scale changes strongly affect both margins. These are fixed-continuation
diagnostics, not a new calibration or proof that timing is generally irrelevant.
The five-page readout is `docs/model/e5f_joint_nested_readout.md` with PDF
`output/pdf/e5f_joint_nested_readout.pdf`. The result README gives source
reconstruction and collection commands. The first-stage cap is closed at
6m35s of allocated-job elapsed time; no additional experiment is running.
