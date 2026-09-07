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


## Full lifecycle and calibration extension (September 6)

The author subsequently requested a full overnight calibration and equilibrium
paths. This extension is isolated on branch `codex/joint-nested-full`. The
first-stage fixed-continuation calculation above is a separate completed
experiment, not evidence that the full model is calibrated.

Every age now integrates the four joint plans before its value is passed to
the preceding age. The original conditional budget and saving solvers,
transaction maps, fecundity and independent child maturation are retained.
The joint plan object is owned by its policy and accepted-price cache.

For a given population, accumulate mass by post-conception family state and
chosen product before applying the original transaction map. Divide this
selected product mass by its post-conception family-state mass to obtain an
effective product kernel. This is an exact computational compression because
subsequent transactions and conditional policies depend on the destination
product and family state. It is valid only for that population. Recompute it
for every date, candidate price and matched birth/control branch; do not use
it as an ex ante household probability or in policy plots. The selected
attempt tenure is preserved for both origin birth/control branches. At the
next date the confirmed-childless control faces the restricted wait menu.

The eleven-dimensional search replaces `kappa_fert` and
`kappa_fert_continuation` by the outer tenure scale kappa in [0.005,10] and a
common dissimilarity lambda in [0.02,1]. Both use log search coordinates. The
upper restriction lambda<=1 is required for the chosen GEV law; the positive
lower bound is a numerical search restriction, not an empirical estimate.
The common lambda across birth orders is explicit and may be rejected by fit
or identification. The prior authorized first-child room-jump upper bound
2.0 is retained. The other nine parameter coordinates and all twelve targets
and weights are unchanged. No target is dropped or demoted.

The complete-loop smoke requires two exact histories (including all target
rows and seventeen PNGs), two small perturbations of every free parameter,
and short versions of all four policy paths. The long controller bounds both
attempted and completed candidates, preserves valid incumbents after declared
price/feasibility rejections, and stops for unexpected source, accounting or
code failures. Final repetitions copy the exact candidate-generator inputs.
The final Jacobian records its own anchor and any missing/one-sided columns;
it cannot certify identification by parameter count alone.

The expectation method remains temporary equilibrium: current prices are
permanent in each household solve and each dated housing market clears. The
post-2023 baseline and supply/LTV/property-tax experiments retain the existing
closed-national finite-horizon closure, with no outside entry, retention one,
2.1 replacement conversion and the exact inherited queue. They are not a
perfect-foresight solution or stationary endpoint. Welfare, rebates, purchase
grants and production promotion are outside this run.

Long-run hard limits and hashes are in the immutable experiment contract.
The target fit, every free parameter/bound, identification diagnostics,
standard graph packets and a short PDF are required deliverables. Safe
numerical completion is not author approval of the nesting structure.

### Inherited-population handoff verification

Each joint-mode dated solution now saves the original input distribution,
before the candidate price applies feasibility projection. The final policy
driver branches every policy from that same original distribution and exactly
replays the fitted feasibility gate as a handoff check. A nonzero projection
is recorded; the original input must not be reconstructed by silently treating
the gated distribution as inherited data. Older experimental checkpoints lack
this field and must be rebuilt. This changes checkpoint metadata only.

### Owner-consumption reporting verification

The unchanged owner optimizer uses actual consumption implied by resources,
housing expenditure and saving. Its legacy output applies a larger reporting
floor after optimization. That can report consumption of 0.040 when a feasible
optimized budget supports only 0.015366, creating an apparent budget violation.
In the experimental joint mode only, feasible solved owner branches now report
\(c=\text{resources}-\text{owner cost}-b'\). The objective, saving choice,
housing product and feasibility requirements are unchanged; infeasible branch
sentinels remain untouched.

The before/after replay at the exact failed supply-policy price verifies fifteen
unchanged value, choice and distribution arrays, unchanged prices, births,
demand and supply, and seventeen identical standard PNGs. Budget-violating mass
falls from \(1.93905\times10^{-7}\) to \(1.16494\times10^{-12}\), below the
unchanged \(2\times10^{-10}\) gate. The standalone checker is
`run_e5f_joint_nested_reporting_check.py`; snapshot e and its receipts are
indexed in the full experiment README. A complete history/policy smoke must
still pass after this source revision. The production solver is unchanged.

### Saving maximization before a full calibration

The subsequent full-loop smoke failed the occupied value check at the 2027
supply-policy state. A same-price, same-population replay exactly reproduces
the local solution; replacing its two conditional saving kernels by the
previously audited exhaustive kernels removes the value decrease. Births
change by -0.03007% at fixed prices. This is a numerical-method comparison,
not a calibrated policy effect or a new equilibrium: the improved policies
require market clearing again. Exhaustive saving must be integrated and
verified through the complete historical/policy loop before launching the
experimental calibration. Its full-history runtime remains unmeasured.
The original gates, economic parameters and empirical targets remain intact.
Read the five-page discussion PDF and receipts indexed in the experiment
README for that earlier diagnostic and review hold.

### September 7: verified saving integration and overnight search design

The author subsequently authorized completion of correctness verification,
followed autonomously by a full experimental calibration and policy paths.
Exhaustive saving is integrated only in the experimental joint-choice mode.
The separate production solver remains unchanged. Two complete starting
histories now reproduce the target fit, parameters, historical path and all
seventeen diagnostic plots exactly; they take about39.4minutes each. The two
all-parameter probes and four two-date policy branches must also pass before
the broad search starts. Their live receipts are indexed in the experiment
README and canonical calibration status.

The reviewed `wide32` controller requests32CPUs and384GB. Its population has
64 parameter vectors: the exact starting vector,47 combinations spanning
eight outer taste scales and six nesting coefficients, and16 nearby joint
perturbations. The remaining nine parameters also vary at every scale-grid
point. Subsequent differential-evolution proposals vary all eleven parameters
against all twelve unchanged targets. The domain and objective are identical
across stages. Sixty-four candidates are stored in two bounded plans, with
global population indices retained across a shared32-worker queue.

The640 attempted-history limit is a ceiling, not a planned completion count:
at the measured starting-vector speed it would require13.1hours even before
policy paths. Each immutable contract records the measured speed, remaining
time and projected capacity. The controller permits up to eight generations
and two local refinement rounds but starts a stage only if all its full
one-hour case-limit waves fit before the search deadline. Four and a half hours are
reserved for22 local sensitivity histories, two exact repetitions and the
four eleven-date policy paths. Total execution is capped at12hours and the
September7 13:35UTC absolute cutoff. Every case writes a heartbeat and saved
diagnostic packet; latest-case and best-so-far summaries remain available.

The wider scheduler accepts the earlier complete smoke only through the
original hash-pinned contract, plans and receipts. It checks the unchanged
scientific source, target system, domain, helpers, numerical gates, closure,
exact historical repetitions and selected-state policy receipt. Independently
executed historical and policy components must be identified as such; they
cannot be relabeled as a completed single controller run. Incomplete policy
finalization cannot produce a full-completion status. A successful search
still does not constitute author adoption of the nesting specification.


**September 7, 03:40 UTC: zero-support diagnosis confirmed; repaired verification running.**
Instrumentation-only job `17090114` completed and reproduced exactly zero
first-birth mass in both stationary comparison branches, with finite valid
probabilities and zero mass difference, at the initial preference guess 0.1062.
This intermediate trial cannot define a conditional birth housing response;
its completed-fertility level remains measurable and can guide normalization.

The isolated repair catches only a typed missing-support exception and records
the auxiliary response as unavailable. Unequal, negative or nonfinite branch
mass still fails. Every normalized-old target row must now explicitly be finite
before the original stationary/dated identity comparison; no row is excluded.
Actual 2019–2023 target support, weights and all original gates are unchanged.
The lead verified the diff and an independent bounded reviewer found no blocker.
Seven operator, twenty accounting, four integration and eleven controller tests
pass. The unchanged old model reproduces all ten reference arrays exactly.

Source `42e0b97` is committed on `codex/joint-nested-full`. The immutable fresh
snapshot `/scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260907e`
has bundle `50b4342797eb271e71b15651c60f4e45c6c740eb6205dd609fcd32501c7428eb`
and contract `9f512dfcaf7dc85b359e0b65d4f8f1543221f4c59f88d82888b91a1ebf4d3e41`.
Failed-case replay `17090361` runs on one CPU/24GB with a 60-minute process cap;
full smoke `17090362` runs on two CPUs/64GB with a 140-minute cap. The latter
requires two exact full histories, two all-eleven-coordinate probes and four
two-date policy paths. At the prior 1,886–2,368-second history and 1,526-second
policy-smoke timings, roughly 100–110 minutes is allowed before queue overhead.
The repaired source has not yet passed these complete checks. Broad calibration
is stopped until it does; the cutoff remains 13:35 UTC, and a new wide contract
must reserve 4.5 hours for final checks/policies. Old smoke cannot certify this
new bundle. A bounded worker is preparing the morning PDF builder concurrently;
the delivered PDF remains the earlier discussion copy until explicitly refreshed.

Evidence: `output/model/e5f_joint_nested_full_20260906a/support_repair_e/`,
`stationary_support_diagnosis/`, and `stationary_support_diff_review.md`.
Job `17089711` remains FAILED, with no new complete calibration case. Production
and the protected manuscript are unchanged. Original c/d snapshots are preserved.


### Small-decline coverage revision prepared, September7 04:40UTC

**September 7, 04:40 UTC: all four fresh histories verified; policy checks running.**
All 84 original artifact hashes pass. Both anchors and both all-parameter probes
exactly match the corresponding snapshot-c target tables, parameters, 253 numeric
historical entries and seventeen standard graphs. Completed histories took
1,970.42, 1,969.84, 1,488.41 and 1,415.85 seconds. The baseline policy passes
2023/2027; three other two-date policies remain in smoke17090362. The small-shock
canary17091265 remains in old-fertility normalization, with a valid bracket now
found. No broad search is running.

A bounded worker and lead check found that the original starting population
omits small preference declines. The prepared controller revision preserves
the exact anchor and47original grid points, replacing16nearby random proposals
with16paired grid proposals having smaller declines. The multiplier is
min(0.5,lambda*kappa/1.6), where1.6 is the anchor's inner taste scale. This is
an initialization heuristic, not an estimated restriction; later search still
varies all11parameters freely. Every target, weight, bound and numerical gate,
as well as the three-timeout stop rule, is unchanged. The paired declines range
from -0.1877535624 to -0.00003064868496. No inference that a target or parameter
region is infeasible follows from this design check.

Experimental source5c70aa1 is committed/pushed. Twelve controller tests pass
locally and onTorch; only initial_population differs in the controller's AST.
The NEW prepared snapshot is `Fertility_Spring26_joint_nested_full_20260907g`,
with68verified files and the same scientific bundle50b434. No contract/search
is created. Use `wide32_support_g/prepare_repaired_wide.py`, indexed in the
experiment README; the old f preparation helper is obsolete. The g helper
also requires a lead canary assessment with verified original evidence before
it will build a contract. Inspect the canary, finish all smoke policies, adopt
or revise the preparation as evidence warrants, then execute full preflight
and launch autonomously. Original e/f source and running jobs are unchanged.
The13:35UTC cutoff and4.5-hour final reserve remain; monitoring is active.


### Complete preflight and broad search, September 7

**September 7, 05:34 UTC: full verification passed; broad calibration submitted.**
Torch job `17093420` is submitted on `cs` for 32 CPUs and 384GB, with a
12-hour Slurm ceiling; the last queue check shows PENDING for resources.
The frozen snapshot is
`/scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260907g`.
Contract SHA is `c5524a0fcf68e78e150da334a98d7b18143af2c9753fa37871a3ada379da36e6`.
Experimental source `7d17162` is pushed; scientific bundle remains `50b434...`.
The 68-file source check, fourteen local/remote controller tests and actual
complete smoke-import preflight all passed before submission.

Fresh smoke `17090362` completed 0:0 in 1h22m10s. Both full-history anchors,
both all-eleven-parameter probes and all four two-date policy paths pass in
one controller. Independent verification covers 84 historical artifact hashes,
all target/parameter/history rows, 170 policy artifacts including all eight
dated pickles, and exact matching of all 360 policy numbers and 136 graphs to
the preceding source. Policy smoke took 1,433.53 seconds, projecting about
2.19 hours for the four full eleven-date paths. No full new policy path or
searched calibration is yet complete.

Small-shock canary `17091265` exhausted its one-hour process cap during old
fertility normalization: fourteen stationary evaluations, best completed
fertility 2.10136049 versus 2.1, above the unchanged 0.0005 tolerance. It has
no normalized-old pass, completed history or valid loss. Its peak 9.253GiB
implies about 296.1GiB for 32 equivalent workers, below the allocation, subject
to monitoring. A timeout does not establish that a target or region is infeasible.

The adopted initialization preserves the anchor and 47 original grid points,
with sixteen paired proposals adding smaller preference declines. Later search
still varies all eleven coordinates freely. The three-consecutive-timeout
threshold still stops new search; active cases now finish within existing
caps, and final verification proceeds from a valid saved incumbent. Unstarted
cases are recorded without a fabricated loss or completion count. Unexpected
scientific failures and rejected required repetitions remain fatal. No model,
target, weight, parameter bound or numerical gate was relaxed.

Search stages must fit before **09:05 UTC**, with **4.5 hours reserved** for
22 Jacobian probes, two exact repetitions and full policies; hard cutoff is
**13:35 UTC**. The contract projects about 192 search histories at the measured
supported-case rate; 640 is an attempt ceiling, not promised completions.
Queue delays, slow cases and the timeout stop can reduce coverage. The monitor
will collect and assess actual results and refresh the morning PDF; the old
delivered discussion PDF remains stale. Evidence, complete smoke fit and
parameter tables, preflight and submission receipt are under
`output/model/e5f_joint_nested_full_20260906a/`, especially `support_repair_e/`
and `wide32_support_g/`. Use task-private `g/tmp` for TMPDIR: the shared
login-node /tmp filled during testing; nothing was deleted. Production and
the author-controlled manuscript remain unchanged.
