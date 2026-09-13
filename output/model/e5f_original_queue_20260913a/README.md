# Original household-queue afternoon experiments

Author approved launch on September 13, 2026. Smoke job **17697809** was
confirmed RUNNING on Torch at approximately 17:08 UTC. No long experiment has
yet passed or been submitted. A successful smoke automatically dispatches the
five independent jobs below. Existing presentation results remain frozen.

## Scientific contract

Reuse the calibrated 2007 structural parameters and raw stationary household
distribution. Use the original household birth-vintage queue throughout:
adjusted births divided by 2.1 become household entrants after four waiting
slots. Seed the stationary adjusted queue from actual age-zero household mass;
retain the separate raw-birth queue. No immigration, observed historical age
reweighting, or subsequent person/headship population-law switch is introduced.

PAYGO payroll tax is 0.179 and pensions balance. Baseline annual property tax is
1%, with equal rebates; historical policy comparisons also use 2%, equally
rebated. Preserve the calibrated supply curve with elasticity 0.63, utility,
parameters, target/weight system and numerical gates.

Standalone IRFs permanently change preference from 0.1489153145785918 to
0.09221854783921073 at once: the entire previously fitted decline. Historical
arms instead refit four successive surprises to the existing four fertility
windows. At each surprise, households expect the new preference to persist;
the final preference is constant after 2023. These are distinct exercises.

All outputs are diagnostics until reviewed. A finite price-path root is not
proof of a stationary endpoint, adequate horizon, existence, or uniqueness.

## Execution and limits

| Stage | Horizon | Per-job numerical budget | Role |
|---|---:|---:|---|
| Smoke | 6 | 40 minutes | Two native no-shock paths, exact root/replay, native diagnostics, original stationary endpoint reproduction |
| Permanent IRF | 24 | 4 hours | Existing finite boundary |
| Permanent IRF | 100 | 5 hours | Existing finite boundary |
| New stationary endpoint + IRF | 100 | 6 hours | Up to 1 hour for terminal equilibrium, then anchored path |
| Fresh historical fit + policies | 6 per forecast | 6 hours | Re-estimate four surprises, then paired rebated taxes |
| Fresh historical fit + policies | 24 per forecast | 6 hours | Same exercise with longer forecasts |

All five experiments run independently after smoke success. Each requests one
CPU and 32 GiB; numerical libraries use one thread. The absolute deadline is
seven hours after preparation, 00:07 UTC September 14 / 20:07 EDT September 13
(exact Unix timestamp in `spec.json`), including
queue time. Node watchdogs enforce this even if the laptop is closed. Each
evaluation writes progress; controller heartbeat interval is one minute.
Incomplete paths and failed roots must retain their diagnostic status.

App monitor `monitor-original-queue-afternoon-experiments` checks every thirty
minutes and reports meaningful changes only. Numerical jobs and conditional
dispatch do not depend on that app monitor or local connectivity.

Root evaluation limits: finite forecasts at most 24 evaluations per root;
new stationary endpoint at most 24 evaluations (16 in smoke); fixed-terminal
100-period path at most 8 evaluations. Historical work uses bounded existing
scalar shock search and forecast controls, with a one-hour policy reserve.
This is not a full structural recalibration. Solve totals are adaptive and are
bounded by these evaluation and wall-clock limits.

## Source and reproducibility

- `code/model/tools/e5f_original_queue_experiment.py`: original-queue adapter.
- `code/model/tools/e5f_original_queue_terminal.py`: stationary renewal/fiscal root and full native one-step audit.
- `code/model/tools/run_e5f_original_queue_experiments.py`: smoke and experiment controller.
- `code/model/tools/build_e5f_stationary_shock_figures.py`: supplementary six-panel IRF figures alongside native diagnostics.
- `code/cluster/prepare_e5f_original_queue_experiments.py`: immutable preparation and smoke-gated dispatch.

Python compilation passed before launch. Native smoke has not yet passed.
The endpoint solver imposes births/2.1 equal to household entry and balances
both fiscal accounts; population scale follows the unchanged housing supply
curve. It does not reuse the person-law endpoint with migration set to zero.

Remote batch:
`/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/afternoon_original_queue_20260913a`

The local `spec.json`, `history_manifest.json`, and `submission.json` are exact
copies of frozen remote receipts. Source and input hashes are pinned there.
`submission.json` records initial smoke submission only; `dispatch.json`, when
created, is authoritative for subsequent job IDs. Read `smoke/summary.json`,
`controller_failure.json`, root receipts and diagnostic metadata before calling
anything successful. Do not overwrite pinned sources while jobs use them.

Collect receipts, JSON/CSV summaries and figures into this folder; avoid large
native checkpoints unless needed for verification. Do not update slides or
replace the existing verified-history packet automatically.
