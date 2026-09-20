# Specification follow-up: September 20–21

## Authority and deliverable

On September 20, Tommaso explicitly authorized an active goal, check-ins,
continued experimentation, and longer cluster jobs while away until tomorrow.
The objective is a reviewed decision packet: a proposed provisional model,
evidence for the few choices that matter, unresolved questions, and the
implementation/measurement/calibration sequence after an author decision.
Experimental results do not adopt a new baseline, target system or policy.

Aim to deliver by **09:00 America/New_York, September 21**. This is a working
delivery target, not a deadline explicitly supplied by the author.

## Stages and initial resource envelope

1. Collect and review the existing Claude Max/Fable assessment. Session
   `fbfab0f8-fd81-4d5e-9453-23b6912efaf9` uses `claude-fable-5-1`, maximum
   effort, read-only tools, 60 turns and a 20-minute wall-clock cap. It began
   at 15:15:34 UTC September 20. Raw receipts are in
   `tmp/fable_specification_review_20260920/` at the repository root.
   A failure does not authorize an automatic restart.
2. Diagnose the credit/lifetime-cohort response using saved outputs first.
   Distinguish changed choices, evolving cohort composition, borrowing-support
   limits and measurement. The identical high-dose outcomes need a separate
   explanation. Do not presume either a bug or an economic sign theorem.
3. Establish which existing housing, ownership, wealth and fertility age
   profiles can be compared with consistent empirical definitions. Historical
   plotting shortcuts are not new calibration evidence. For example,
   `house_size_age_model_vs_data.py` uses a fixed renter-room proxy and cannot
   provide an exact realized-room comparison without correction.
4. Use these findings and Fable's assessment to specify at most two additional
   model variants worth testing. Changes must isolate a stated economic
   question and have an unchanged-source control. Limited conditional fitting
   can test feasibility; it is not serious calibration or adoption.

Initial envelope: **at most 16 new fixed-price lifecycle/cohort evaluations
and 32 new full stationary evaluations**, including controls, loop smokes and
selected-point repetitions. Replaying saved arrays does not count as a new
solve. At most eight single-threaded workers concurrently. Each batch must
state its actual dimensions, solve count and timing estimate from the latest
relevant measurement; these maxima are not an instruction to fill the budget.
No new batch after 02:00 September 21; planned numerical work must finish by
07:30 to leave time for collection and review. Each submitted allocation and
controller must enforce its own smaller applicable cap. Do not expand this
envelope silently.

## Required launch contract

Before submission, record a hypothesis, exact source/input/target hashes,
parameter overrides, population/entry contract, held-fixed objects, outputs,
per-case and total time budgets, case limits, and stop conditions in a launch
manifest in this directory. Estimate time from an observed comparable solve.
Use a fresh immutable remote snapshot and the exact-loop smoke; production
depends on successful smoke through Slurm, without laptop-dependent chaining.
No production job is authorized merely by appearing in a draft plan.

Keep the original objective, target definitions, weights, bounds and numerical
gates pinned for any comparable fit. A model variant may require a different
measurement contract; if so label it explicitly and do not compare losses as
though identical. Never add a free parameter without identifying variation or
an external restriction. All core numerical changes require lead review
against the model equations before they are used.

Progress/checkpoint output is required every case or five minutes, with a
latest-completed summary and best-so-far summary when selection is involved.
Investigate after 30 minutes without a heartbeat/checkpoint. On failure,
preserve evidence, stop the affected stage and cancel only never-started
dependents. A second attempt requires a documented changed hypothesis or
method; do not relax gates, silently change populations, or repeat a failed
search with a new label. Independent healthy branches may continue.

Each solved case must retain the standard 17 diagnostic plots and numerical
checks. A selected fitted candidate needs the complete target/parameter tables,
actual search bounds, source/objective fingerprints and two exact repetitions.
Keep snapshot total/first birth flows separate from explicit lifetime cohort
births, lifetime first-birth probabilities and the model's stationary
normalization. Across-checkpoint comparisons do not hold prices, preferences
or entrant wealth fixed unless explicitly verified.

## Check-ins and present state

### 12:34 Eastern: credit replay reproduced; income-grid jobs submitted

Corrected credit replay **18078912/18078913** completed all three cases. The
lead checked every old-scalar/cohort reproduction gate (1e-10), source checks,
zero occupied value drops, zero excess-budget mass, and 17 plots per case.
Case runtimes were 61.2, 55.1 and 54.4 seconds. The two larger credit allowances
have bitwise-identical snapshot populations and birth flows, but their complete
policy arrays differ. Thus identical outcomes do **not** mean identical policies
over the full state space. Occupied-state comparisons are being collected before
interpreting this difference. [Raw policy comparison](credit_policy_retention_v2/results/policy_comparison.json).

Income-grid smoke **18079046**, dependent production **18079047**, are submitted
under [income_grid_cohort_v2/submission.json](income_grid_cohort_v2/submission.json).
Frozen source and the nine runtime helpers are copied and checked before adding
the experiment driver; every staged Python file is pinned, inputs are immutable,
and the production driver checks the completed baseline policy/cohort/plot gates.
Allocated household evaluations so far: four credit replays including the failed
v1 receipt attempt, plus three income-grid cases, seven of sixteen. No stationary
equilibrium evaluation or rental-wedge solve has yet been allocated.

### 12:30 Eastern: replay failure isolated; corrected chain and grid staging

Credit replay v1 smoke **18078707** failed after its one household solve,
cohort and all 17 plots completed, during final import-origin receipt assembly.
The check rejected Python's `__mp_main__` alias even though it pointed to the
exact pinned driver. Production **18078708** was cancelled by dependency.
The failed attempt and artifacts remain in `credit_policy_retention_v1/`.
The changed method accepts that alias only at the exact pinned driver path;
a focused fixture verifies that foreign paths are still rejected. No economic,
source-hash, population or numerical gate changed. Fresh snapshot
`credit_policy_retention_v2` has smoke **18078912**, production **18078913**.
Including the completed v1 smoke, this replay chain uses at most four household
solves. Its full arrays remain on Torch for compact comparative collection.

The original `income_grid_cohort_v1/launch_manifest.json` is an unsubmitted
implementation draft. The lead replaced its incomplete staging logic with
explicit frozen-source copies and pinned runtime overlays. The actual launch
uses fresh `income_grid_cohort_v2`; submission status is its `submission.json`
when present. Design remains 5-by-3 smoke, then 9-by-3 and 15-by-3; three solves,
900 seconds each, one CPU/24 GB, 20/40-minute allocations, 1050/2250-second
stage process limits. Larger grids use the same three-point transitory process
and native conditional-entry rule; their entry wealth marginal may change.
No grid change is adopted, and lifetime-cohort results are not stationary fit.

The isolated wedge port now has a segment-by-segment independent saving oracle
and is being made reproducible as a frozen-source patch. Claude Max is doing a
bounded two-file driver/launcher pass; lead review and full controls still gate
launch. No wedge household solve has been launched.

### 12:17 Eastern: saved-array results reviewed; credit replay submitted

Housing collection **18078558** and saved-cohort debt jobs **18078581/18078582**
completed with all retained gates passing. Housing outputs reproduce checkpoint
uncapped rooms and ownership exactly, then apply the empirical cap at nine before
weighting. The reproducible overlay is in
[housing_profiles_v1/reproducible_overlay_v1](housing_profiles_v1/reproducible_overlay_v1/overlay_receipt.md).
The model room profiles remain high into old age, where the empirical profile
falls; ownership is too low earlier in adulthood. These are descriptive
comparisons across checkpoints with different prices, preferences and entry
wealth, not causal income comparisons. The model has no DUE structure category.
Complete calibration tables remain in the [overnight readout](../overnight/final_search/readout.md).

The [saved-credit comparison](saved_credit_v1/results/comparison.md) covers all
12 family/dose cases. In the refit, the saved cohort populations for credit doses
1 and 5 are bitwise identical; in the pilot they differ by less than 1e-46 in L1.
Neither has observed cohort mass at the lower wealth-grid boundary. This weakens
the particular grid-floor explanation, but retained arrays cannot establish
native borrowing-constraint binding or policy identity. Negative liquid wealth
is not itself a measure of unsecured borrowing. Shares aggregated across ages
weight each age by surviving cohort mass.

Credit policy retention smoke **18078707** and dependent production **18078708**
are now submitted. [Submission](credit_policy_retention_v1/submission.json) and
[immutable manifest](credit_policy_retention_v1/launch_manifest.json) govern the
fresh remote `credit_policy_retention_v1` snapshot. One baseline solve precedes
two treatment solves, with 10/20-minute allocation caps, one CPU/24 GB and
per-case 600-second outer limits. All original source, population, numerical,
cohort and 17-plot gates remain. The unchanged replay must reproduce the old
case scalars within 1e-10 and retain full policy arrays. The comparison of the
policies is an experimental result, not a gate requiring them to be identical.
Three of the sixteen household evaluations are allocated to this chain.

The deterministic [income-grid resolution diagnostic](income_grid_resolution_v1/README.md)
is complete: ten grids, all stationary/normalization/log-moment gates, baseline
fingerprint reproduced, and two supplemental plots. Its discretization limit
is the continuous endpoint process. Its proximity to the exact annual block
variance at a particular finite grid does not validate the aggregation mapping.
The separate three-solve cohort sensitivity launcher and isolated rental-wedge
port remain under review and have not been submitted.

Two earlier housing collection attempts failed before model-data collection:
18078441 used the wrong Python environment; 18078503 staged under login-only
`/tmp`. The successful method uses explicit Anaconda and shared `/scratch`.
No household solves were consumed by these collection attempts. The credit
launcher was repaired through a bounded Claude Max review before lead review;
no active model core or manuscript was edited.

### 11:52 Eastern: launched earnings diagnostic; housing source replay passed

Torch jobs **18078286** (smoke) and **18078287** (full, `afterok` smoke)
use `income_aggregation_v1/submission.json` and `launch_manifest.json` here.
Remote root is
`/scratch/td2248/projects/Fertility_Spring26_specification_20260920/income_aggregation_v1`.
Results are under the corresponding repository-relative
`output/model/native_financing_diagnostic_20260919/specification_followup/income_aggregation_v1/results/{smoke,full}`
inside that root. The full design is 20 independent batches of 20,000 annual
income paths, 120 years each; smoke uses two batches of 2,000 paths, 40 years.
Seed is 20260920. Internal caps are 60/720 seconds, process caps 540/840
seconds, allocation caps 10/15 minutes, one CPU/4 GB each. Scientific gates
test mean normalization and analytically known level covariances; neither a
successful run nor small approximation gaps adopt a new income process.
The household/stationary solve budgets consumed by this chain are both zero.

`housing_profiles_v1/full/target_recomputed.json` reproduces all four active
housing targets exactly after the one-chunk smoke. The full empirical pass
traversed 5,848,121 raw records in 24 chunks and 14.24 seconds. No model solve
or target change occurred. National figures must retain the same row-specific
definitions: ownership ages 30–55 in DUE structures is 0.6762604; all-age
ownership is a different statistic. Supplemental empirical age plots and
saved-array model comparisons are the next measurement step.

The proposed wedge experiment is still blocked on correct rental expenditure
in budget audits and isolated-source reproduction. Existing active switch-off
tests are insufficient to establish reproduction of the frozen checkpoint.
The saved-credit diagnostic is being narrowed to actual cohort debt/support
and policy identity, avoiding unsupported native constraint-binding claims.

### 11:54 Eastern: earnings diagnostic completed and reviewed

Both jobs completed with exit zero (7/8 seconds). The full result passes every
normalization, covariance, batch-completion, source-hash and plot gate. Results
and the collection receipt are in `income_aggregation_v1/`. The lead inspected
the covariance figure and the full numerical receipt.

The current 15-state chain matches the continuous endpoint approximation's
log covariances to numerical precision, but not its level covariances. At lag
zero, level variance is 0.961758 (15-state), 1.204895 (continuous endpoint),
and 1.155236 (exact four-year average of annual earnings). At lag one the
corresponding covariances are 0.673591, 0.848034 and 0.850119. Monte Carlo
reproduces the analytically known block moments within its declared uncertainty.
This establishes a numerical-distribution difference, not its effect on
household choices, calibration or the preferred specification.

A bounded deterministic grid-resolution diagnostic now compares persistent
node counts 5/7/9/15/25 and transitory counts 3/5, holding annual parameters and
period mapping fixed. It must reproduce the existing 5-by-3 payload and report
both log and level moments. This uses no household solves or new Monte Carlo.
A household sensitivity run remains conditional on reviewing a coherent
population/entry mapping and its actual solve count.

### Reviewed Fable assessment and current work

The initial Max/Fable review completed in 517.5 seconds. A separately bounded,
focused clarification completed in 208.1 seconds in the same session. The
[assessment and clarification](../../../../docs/model/structural_model_specification_fable_20260920.md)
are retained verbatim, with the [review receipt](fable_review_receipt.json).
The clarification withdraws several initial overclaims. Counted moments do not
establish identification; bound hits do not prove target incompatibility; the
no-type candidate cannot be evidence that its nonexistent types caused the
wealth dispersion; a finite rental wedge only approaches a hard cap in a limit.
The plateau doses are four and twenty annual earnings amounts. Estate
liquidation and deterministic tenure remain substantive choices.

The lead does **not** accept importing BGM's transitory-variance correction into
the project's different income concept, or arbitrary economic success cutoffs
from the clarification. These are recorded suggestions, not experimental facts
or adopted restrictions. The next candidate intervention is the rental-size
price slope, with the owner premium retained and intercept zero; a cap-10,
zero-slope control must separate added rental support from price effects.
Implementation and launch remain conditional on source/accounting review.

Three bounded deliverables are being implemented by Luna with separate file
ownership: saved-array credit/debt/support diagnosis; exact-moment and simulated
annual-to-four-year earnings aggregation diagnosis; and early-ACS housing age
profiles. They do not themselves solve a new household model or change targets.
The earnings exercise measures the existing aggregation approximation and does
not silently install a new process. The empirical exercise uses the actual
active **2005–2006, 42-metro, cap-at-9** source, and separately labels national
diagnostic figures. The older 2012–2023 cache is metro-restricted and is not an
appropriate national or active-period replacement.

Current target provenance is the early-housing builder and candidate packet in
`output/model/e5f_matched_pf_20260909a/design_research/housing/`, with the exact
working target contract in `initial_calibration_contract/working_weights.csv`
under that experiment. The historical August 17 room receipt does not contain
the current capped-room values and must not be substituted.

The existing `earnings-refit-follow-up` heartbeat has been repurposed for this
goal. It checks every 30 minutes, gives a concise author-facing check-in about
every three hours between 08:00 and 22:00 Eastern, and reports substantive
findings, failures or needed decisions promptly. Otherwise unchanged polling
stays quiet. Last author-facing start update: approximately 11:22 Eastern,
September 20. End the follow-up after the final reviewed packet is delivered
and every submitted job is terminal or explicitly blocked.

Historical startup: at 11:24 Eastern, Torch authentication succeeded and the user queue was empty.
Three bounded Luna tasks cover credit-code grounding, empirical-profile
inventory, and deterministic collection of Fable's response. At that time no new model
job had been submitted. The completed September 19 overnight jobs are
historical evidence; their failed/cancelled jobs must not be revived.

Keep task-owned changes separate from the repository's substantial unrelated
dirty work. The active decision ledger and author manuscript are not part of
this editing scope. Update this file and `CALIBRATION_STATUS.md` with reviewed
launches/results; commit and push only task-owned changes.

## Numerical follow-up design (submission status above)

The source and array inventory changes the next steps. Factorial arms retained
cohort distributions but **not treatment policy arrays**. Baseline policy and
current-population arrays remain available in each pinned checkpoint. Thus the
credit saved-array diagnostic now studies cohort debt/support and records policy
identity as unavailable; the housing collector reads baseline checkpoints.
No missing arrays are silently reconstructed.

Two small fixed-price job chains are being prepared for lead review:

- **Credit policy retention:** refit baseline control smoke, then credit doses
  1 and 5 at financed share 0.8 and rental cap 6. Exactly three full household
  solves, using the unchanged numerical path, but now retaining treatment policy
  arrays. Every result must reproduce its previous arm's scalar/cohort receipt.
  This is a disclosed replay with new retained outputs, not a restarted search.
- **Income grid cohort sensitivity:** refit 5-by-3 baseline control, then 9-by-3
  and 15-by-3 income grids, with the same annual parameters and period mapping.
  Exactly three full household solves. Price and all non-income parameters are
  fixed. Native conditional entry is rebuilt at each grid, and its wealth/income
  marginals must be reported. Lifetime-cohort statistics are not stationary
  calibration moments or direct ACS fits. A larger grid is not adopted by this
  exercise. The apparent match of the 25-by-5 grid's level variance to the exact
  annual block variance reflects opposing approximation differences; the
  discretization convergence target is the continuous endpoint process.

Both chains need immutable source/input manifests, baseline reproduction,
full numerical gates, standard 17 plots, and Slurm smoke dependencies before
submission. Proposed caps: credit smoke 10 minutes/production 20 minutes;
income smoke 15 minutes/production 60 minutes. Their combined proposed six
household solves leave ten of the initial sixteen available. The rental wedge
port is still under development and no wedge solve is yet counted/submitted.
