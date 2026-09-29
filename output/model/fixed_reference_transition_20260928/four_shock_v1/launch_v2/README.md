# Bounded shock estimation launch

Tommaso authorized launch after readiness on September 28, and read-only overnight
monitoring for major failures. This supersedes the earlier preparation-only limit.
Cheaper workers own implementation; the lead reviews economics, code and evidence.

**Submitted:** four successive surprises **18765931**, one permanent shock
**18765932**. Both require successful smoke job **18765176** (`afterok`), with
invalid dependencies cancelled by Slurm. They were queued, not yet fitting, at
submission. `launch_receipt.json` records source pins, budgets and job IDs.
**Final state:** smoke **18765176 PASS**, 43:58 wall time. Both fit jobs then
failed within seven seconds at the driver's 2-GiB cache guard, before estimation.
The launcher requested 64 GiB. The monitor detected this and was paused; no
repair or restart occurred under the overnight monitoring-only instruction.

The reference is the immutable block0506 September 28 export. The only intended
economic change is the author-requested fertility-preference shock: four successive
permanent surprises, or one permanent 2007 shock fitted to the 2020–2023 window.
Earnings, initial wealth/income, credit, inherited distributions, both entry queues,
payroll tax and fixed physical housing stock retain the prepared contract. The
inherited estate settlement remains provisional. This is not a policy simulation
or adoption of another chat's credit/two-birth experiments.

## Readiness

Torch job **18765176** runs the reviewed source under `source/code/`, with one CPU
and 24 GiB. It runs the pure tests, measures five baseline residual mappings at
12 dates, solves the unchanged six-date equilibrium, replays its first date, and
checks the remaining five dates from the carried distribution and both queues.
The independently pinned unchanged endpoint is reused from job 18759951.
No historical fertility is fitted by this smoke test.

The conservative ceiling is 156 Bellman calls before exact caching. Previous
unchanged six-date mappings took about 73 seconds with one solve and eleven cache
hits. Perturbed derivative mappings require additional policies: allow 50 minutes
for all five seed mappings and 60 minutes total, with a 70-minute Slurm ceiling.
This estimate is uncertain; the fixed call count and deadlines are the stop rules.

`smoke/readiness.json` must report PASS with the current nine source pins. The
measured Jacobian can then be reused by both fits with pinned provenance, avoiding
duplicate derivative measurement. Native nonlinear residuals, terminal checks and
fresh replay determine acceptance; the approximate Jacobian only initializes roots.

## Submitted estimation budget

Both cases will compare 104- and 128-date forecasts. Each scalar candidate solves
its own stationary endpoint, then both conditional perfect-foresight paths. There
are at most 12 candidates per fitted shock, 24 endpoint evaluations and 16 mappings
per horizon. This implies conservative ceilings of 357,728 Bellman calls for the
four-shock case and 89,528 for the one-shock case, before caching and time stops.

Historical 104-date mappings took roughly 29–57 minutes in an older model; those
are planning evidence, not a speed measurement for the current shocked model.
The full evaluation ceiling would therefore exceed an overnight allocation.
Each fit has a 12-hour total budget, ten hours per candidate, six per path,
one per endpoint, and 90 minutes per mapping; Slurm stops it after 13 hours.
Partial progress is saved; an interrupted or unconverged fit is not an estimate.

Each job requests four CPUs and 96 GiB, with one estimator process and one
numerical thread. Torch rejected the one-CPU/96-GiB `cs` request before creating
a job; its test-only scheduler check accepted four CPUs with that memory.
The 64-GiB exact cache can hold a 128-date collection at the observed baseline
policy size (about 364 MB each), avoiding eviction before the forward replay.
Cache hits and actual solves remain recorded, so this capacity calculation is
not presented as a measured shocked-path speedup.

The original 46 tests passed inside the smoke. The updated seven-test controller
suite also passed on Torch, including its explicit 64-GiB plan setting. Shell
syntax and fail-closed rejection of an invalid controller pin were verified.
`source/` retains the smoke sources; `source_launch/` preserves those nine model
sources and adds the reviewed cache-budget option and gated Slurm launcher.

Only launch after the smoke passes. Overnight monitoring reads Slurm state and
small receipts, and reports major failures. It must not edit code, loosen gates,
change parameters, cancel, restart, resubmit or launch additional experiments.

## September 29 bounded audit — no refactor or relaunch

One cheaper worker reviewed the launcher interfaces; the lead verified findings
and the saved plans on Torch. No model was loaded or solved, and no source changed.

- **Required launch correction:** replace the inconsistent cache setting/limit
  with one resource contract. `validate_plan(..., launching=True)` rejects both
  actual plans at estimator line 135. Changing only the cache to 2 GiB in an
  in-memory diagnostic copy passes the remaining structural checks; this was
  neither a saved plan edit nor a recommendation to sacrifice cache capacity.
- **Required preflight correction:** controller construction and its cache test
  do not invoke the real validator. The driver's inspection command at lines
  562–563 also skips execution checks: both bad plans report no missing settings.
  Use the same full, no-solve validation before submission and at execution.
- **Required retry preparation:** gated launcher lines 33–42 refuse existing
  `plans/<kind>.json`. These files now exist. A retry needs a fresh run location;
  preserve overwrite protection and the failed evidence. Source edits also
  require appropriately updated readiness/seed provenance, not copied claims.
- **Small performance correction:** `fit_one` samples the initial level at line
  50 and again through the root's unconditional initial sample (root line 128).
  Reuse the already certified initial observation; retain the independent final
  replay. Otherwise each stage repeats an entire candidate's two forecast roots.
- **Budget clarity:** endpoint/path deadlines are tested between evaluations;
  `guarded` only includes the mapping, candidate and total deadlines. An ongoing
  call can overrun its stage limit. Pass the stage deadline into the existing
  watchdog if those limits are intended to be hard.

Both current source fingerprints, readiness receipts, complete target contracts,
and derivative receipt/matrix fingerprints agree. The shell environment's
identification-literal check is redundant, but actual contract authentication
exists in `load_reference` lines 221–232; it is not missing. Keep market,
population, timing and genuine source-identity checks. The audit finds no further
current-plan startup mismatch after isolating the cache cap; it does not certify
long shocked-path convergence or constitute a whole-codebase review.
