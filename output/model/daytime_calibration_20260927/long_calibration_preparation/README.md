# Long calibration preparation — September 27

**Preparation only. No job, model solve, runnable contract, or specification change.**
This is a proposed eight-hour Torch design for the lead to finalize after the
frictionless-transition checks and remaining author decisions. Eight hours is a
planning option, not an imposed deadline or a new launch authorization.

## Starting evidence

The best independently repeated candidate has loss **41.112**, from
`../local_continuation/run_v1/selected_export/`. The original overnight benchmark
has loss **42.282**, from the supervised overnight main continuation, `de_0093`.
Keep both: the new candidate is a better starting point, while the ongoing credit
benchmark/transition uses the original benchmark. Do not silently replace its
initial distribution or preferences with the new candidate.

The existing contract has fourteen displayed restrictions: thirteen scored
moments and completed-fertility normalization. Nine coordinates are searched;
the normalized child-benefit level is the tenth fitted parameter. Current source
inventory hash is `49566b92deae0d27014ebe4204cc88209ea27a664e90ca697d2e377d2f8f4737`;
current target-and-weight fingerprint is
`e79b1d277510e21f7fdf6ac656aca37cd097d9f523500de7dff9ca6f11132f27`.
These identify the current comparison system; they are not permission to keep
it after an author decision changes the specification.

Local objectives took roughly five to six minutes, typically six stationary
solves. These are Mac measurements, not a certified 24-worker Torch benchmark.
Two fresh exact-loop smokes and a full-width first batch must establish cluster
speed/memory and cross-host agreement before the main search proceeds.

## Proposed stages and ceilings

All stages share one fixed global deadline. Queue delay, smokes and preparation
must be distinguished from numerical runtime. Pin the numerical start/end once
an allocation exists; never extend a running deadline. One allocation can avoid
inter-stage queue delays. Twenty-four single-thread workers are the ceiling,
not a promise that the scheduler will provide them.

| Stage | Window from numerical start | Objective ceiling | Purpose |
|---|---|---:|---|
| Acceptance | 0–30 min | 4 | Two exact-loop repeats at each retained starting candidate; authenticate full tables/checkpoints and plots. |
| Broad starts | 30–90 min | 96 | Diversify starting points in explicitly approved transformed coordinates; include retained anchors. |
| Parallel searches | 90–240 min | 576 | Three independent eight-worker populations/streams, up to 24 rounds each, preserving genuinely distinct starting regions. |
| Local refinement | 240–300 min | 192 | Refine the best separated candidates; include half-step checks for nonlinear coordinates before any Jacobian-based move. |
| Sensitivity and profiles | 300–390 min | 126 | Up to 54 finite-difference evaluations and 72 conditional-profile evaluations, separately labeled from search. |
| Verification/export | 390–480 min | 6 | Freeze selection; two independent repetitions for each of at most three finalists, full standard plots and report tables. |

Maximum: **1,000 full objectives**, approximately **6,000 stationary solves** at
the observed six-solve count; the actual normalization algorithm permits up to
23 solves per objective, so 6,000 is an expectation, not a hard solver-call bound.
At six minutes and 24 continuously busy workers, 1,000 evaluations require about
250 minutes; generation barriers, tails, acceptance and reporting use the rest.
The stage and global deadlines bind before the case ceilings. Do not launch
extra cases merely because workers are idle. A slow initial cluster batch should
reduce remaining case counts within the same deadline, not dilute tolerances.

These ceilings require a new reviewed orchestration layer. The existing frozen
controller permits at most 30 rounds of 24 cases (720 search cases) per contract;
do not pretend that it implements this entire staged plan unchanged.

## Code to reuse, and what must be implemented

- `code/model/tools/run_e5f_utility_overnight_calibration.py`: authenticated
  worker evaluation, complete-target validation, candidate identity, numerical
  failure classification, exact-repeat comparison and verified 17-plot export.
  Its existing search is random coordinate/joint perturbation, **not differential
  evolution**, regardless of `de_` case labels.
- Frozen `tools_v4/run_e5f_utility_comparison_search.py`: owned process groups,
  hard child deadlines, finite batches, immediate fatal-stop dispatch rule and
  preserved censored/inadmissible results.
- `code/model/tools/run_e5f_local_continuation.py`: useful small-plan pattern for
  explicit point lists, pinned anchor files, five-second heartbeat and exact
  anchor checks. Its two-worker/29-minute limits must not be bypassed for Torch.
- `code/cluster/run_e5f_supervised_lane.sh` and the daytime preparation/acceptance
  tools: reuse smoke-then-acceptance sequencing and source/table checks. The
  September 27 daytime submission was cancelled before allocation; no result
  from that job certifies 24-worker operation.
- `run_e5f_simple_fertility_overnight.py` contains an older genuine DE design,
  but hard-codes an obsolete 11-coordinate/12-row system and historical cases.
  It is a design reference only; never submit it for the current model.

Implement the staged candidate generator/controller in a fresh immutable
snapshot, with synthetic dispatch/deadline/fatal/restart tests and exact-loop
smokes. Choose explicitly between independent random-search streams and genuine
population-based DE; the latter requires tested mutation/selection semantics.
Checkpoint RNG state, population/incumbents, proposed/completed case identities,
remaining evaluation budgets, and all absolute deadlines. Restarting requires a
new linked controller receipt with only unconsumed time/cases; no silent reruns.

## Decisions required before a runnable contract

1. **Baseline credit rule:** current ongoing collateral constraint or an
   author-adopted DUE-style treatment of inherited owner debt. DUE is not adopted.
   Payment-to-income remains deferred. The frictionless benchmark is a diagnostic,
   not the default economy to calibrate.
2. **Other fixed inputs:** close or explicitly retain/defer the 2% interest rate,
   tenure-choice scale, rental menu, income/entry inputs and their measurements.
   Finer housing grids and conception schedules remain separate unless adopted.
3. **Objective and identification:** explicitly retain or change the provisional
   early-fertility weight/curvature interval and measurement approximations. Do
   not drop targets, silently substitute identity weights, or pool losses across
   incompatible systems. Alternative weights require separate lanes and common-
   weight comparisons, with named identifying moments for affected parameters.
4. **Source choice:** frozen calibration source versus the newly integrated native
   production source. If using native production, repeat the exact current best
   candidate with all fourteen target and 31 parameter rows; existing exact
   native replays at the original benchmark do not certify the improved point.
   Hash the complete source/input/runtime bundle and accepted economic changes.
5. **Search design:** approve duration, coordinates, initial-region coverage,
   proposal scales, per-objective cap and actual algorithm. Avoid uniform raw
   sampling over the entire very wide taste-scale bounds without a reasoned
   design. Wider exploration can retain bounds without changing them.

The first Jacobian is local and step-sensitive: beta, H0 and chi show substantial
nonlinearity. It does not justify an unrestricted Newton update or certify
identification. Half-step diagnostics should precede their gradient-based use.
A conditional profile must state which other parameters are reoptimized; a
one-at-a-time slice is not a profile likelihood or a confidence interval.

## Monitoring and stop rules

Keep latest-completed and best-so-far summaries, every complete 14-row fit and
10-fitted-parameter table, raw 31-parameter tables, case logs and heartbeat.
Track actual worker count, queued/running/completed/inadmissible/timeout/fatal
cases, checkpoint age, bounds and economic misses. At least hourly inspect the
unchanged seventeen plots for a saved candidate, plus separately labeled
occupied-region policy/distribution diagnostics. Flag the known ownership
nonmonotonicity rather than accepting a low loss as policy certification.

No progress for thirty minutes triggers bounded diagnosis. Contract/source,
accounting or unexpected failures stop dispatch immediately; demonstrated fixes
need a new versioned snapshot and tests, leaving old evidence untouched.
Proposed search stop: two consecutive batches with at least half inadmissible
or three all-timeout batches; retain valid incumbents and investigate. All
stages obey finite ceilings; stagnation prompts a recorded stage decision,
not an automatic deadline extension. Reserve final repetitions/export even
when search is disappointing. Never claim convergence solely because a budget
expired or a repeated candidate is reproducible.
