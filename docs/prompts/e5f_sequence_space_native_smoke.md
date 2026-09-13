# Codex worker task

## Goal

Connect the existing sequence-space prototype to the actual frozen Torch model
and attempt one tiny native check. The author approved a 15-minute Terra task,
including permission to proceed despite the previously requested weekly usage
floor. Stop and report at the stated deadline; do not broaden or restart.

## Scope

Exclusive ownership: `code/model/tools/e5f_sequence_space_prototype.py`, its
`test_e5f_sequence_space_prototype.py`, and `docs/model/e5f_sequence_space_prototype.md`.
If needed, add only `code/cluster/check_e5f_sequence_space_native.py`.
Results: `output/model/e5f_sequence_space_prototype_20260913/native_smoke/`.
Write a compact `worker_report.md` there, including any job ID and exact remote
output path. No edits to the active scientific kernels or the other jobs.

## Context

Use full startup with bounded leading excerpts from mandatory context. The
first worker created an interface and four toy tests, but no native mapping or
official-package factorization ran. Its `--native-smoke` only prints a command.
Implement actual native evaluation rather than treating that stub as a test.

The native root receives house prices in levels. The prototype explicitly
uses log prices, exponentiating at its native boundary; pensions and rebates
remain in levels. Keep that coordinate transformation explicit. The native
residual is housing relative imbalance and 200 times each of PAYGO and rebate
relative imbalance, stacked by equation type. Never infer fiscal residuals
from unscaled currency amounts. Read the authoritative `dated_residual` and
`stack_dated_residuals` in frozen `e5f_rebated_surprises.py`.

Torch spec:
`/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/afternoon_original_queue_20260913a/spec.json`.
Import its `source/run_e5f_original_queue_experiments.py`, call
`load_context(spec)`, and set `c.spec_path`. This verifies scientific hashes.
Read `code/cluster/replay_e5f_original_queue_diagnostics.py` for a working
connection to that frozen runtime, original queue, cache, household audits,
fiscal records, and pinned inputs. The native path has T+1 values including the
terminal value; it has T dated quantity rows. Put Slurm logs outside the output
directory if the worker refuses a nonempty output directory.

Use a two-period no-shock path around the original 2007 stationary equilibrium:
native inherited households and both four-vintage queues, original parameters,
constant original price/pension/rebate, and original stationary continuation.
No immigration, no age reweighting, no population normalization during replay.
PAYGO tax .179, equal property-tax rebate, existing supply curve and all gates.

## Do not touch

No production code, calibration, target weights, numerical gates, parameters,
population law, running sources/jobs, slides, manuscript, commits or pushes.
No new large run or dense Jacobian. No global environment/package changes.
Do not claim a fast-news Jacobian, equilibrium convergence, or speedup from a
finite-difference directional check. This is a connection and measurement test.

## Required output

A working native evaluator, a small recorded test if possible, and a concise
report separating native evidence, official-package evidence, and missing
derivative construction. A matrix factorization on a toy system is package
API validation only. If the native job is still running at the worker deadline,
return its ID and output directory immediately; the lead will collect it.

## Verification

At most ONE Torch job, horizon 2, 1 CPU, 24 GiB, 10-minute Slurm cap, numerical
threads one, at most SIX full native mappings. Account `torch_pr_570_general`;
Python `/share/apps/anaconda3/2025.06/bin/python`. Use an isolated, new batch
with immutable source/input hashes. Never run the model on the login node or
Intel Mac. Reuse the already-passed native loop prerequisites.

Suggested six mappings: baseline, exact baseline replay, plus/minus one small
direction at step h and h/2. Use the original scaled residual equations and
native household/mass/feasibility gates. Baseline must meet 2e-4; exact replay
must agree to 2e-10. Report derivative step sensitivity without relaxing it
away or implying global differentiability across discrete choice thresholds.
Keep every native state and queue object in the bridge.

Check the official `sequence-jacobian` package in an isolated target. Inspect
its actual constructor/update interface; do not assume the untested prototype
call is correct. Spend at most two minutes on package acquisition. If it is
unavailable, still connect and test the native evaluator and clearly report
that official-package integration remains unverified. The prior worker's
restricted network failure is not evidence that Torch cannot be reached from
this in-thread worker's tools.

## Stop and report if

The 15-minute deadline is reached, the native scientific contract cannot be
preserved, or a broader rewrite is required. Return concrete artifacts and
the unresolved step. The lead will review every numerical change against the
native equations before trusting any claimed result.
