# Original household-queue afternoon experiments

**20:24 UTC — corrected graph collected and visually verified.**
Replay17702691 exactly reproduced all100 iteration3 rows. The requested
native fertility-rate graph, including the pre-shock steady state, is at
`terminal_restart_v1/fertility_replay_iter3/output/irf_fertility.png` and `.pdf`.
Fertility drops2.10 to1.68186 on impact and rises to2.08277 by the final
period. The carried year400 household mass0.387755 is10.96% above the
stationary0.349454. This is an unconverged iterate. The main run has moved to
iteration5, housing gap1.70615%; its later progress is separate from this plot.

**19:44 UTC — fertility figure:** diagnostic replay job17702691 is running
on cs609. It reproduces saved 100-period iteration3 with identical controls,
saves native fertility rates and terminal state gaps, and automatically plots
the fertility-rate panel with the initial steady state. Local contracts and
collection instructions are in `terminal_restart_v1/fertility_replay_iter3/`;
remote directory is the main batch plus `_fertility_replay_iter3`. The thread
monitor includes it. Original equilibrium jobs continue independently.

**19:17 UTC collection:** ten-period job17700926 finished at its eight-mapping
budget without convergence. `terminal_10/` holds its native rows, fertility,
stationary references, root receipt, terminal distances and visually inspected
supplemental plot (regenerate with `code/model/tools/build_e5f_fixed_terminal_horizon_figures.py --case-dir output/model/e5f_original_queue_20260913a/terminal_10`).
After40years household mass0.917922 remains162.67% above terminal0.349454.
The plot labels the failed root explicitly. The stationary100 arm has completed
three mappings, latest score8.48344 versus2e-4 required, and is still running.
Appendix algorithms are in `latex/appendix_solution_algorithms.tex` with a
standalone two-page extract at `output/pdf/solution_algorithms.pdf`.

**18:41 UTC — new author-requested short test:** job17700926 reuses the verified
new steady-state endpoint for a ten-period (40-year) transition. It first runs
the same ten-period root/replay with no shock; only a passed no-shock check
permits the shocked stage. Eight mappings per root, ten-minute smoke ceiling,
45-minute total numerical budget,50-minute Slurm cap, original global deadline
retained. The population is carried forward and its endpoint discrepancy is
reported, never removed by rescaling. Native diagnostic measurement of period
fertility is retained; aggregate birth counts are not substituted for a rate.
Additional plot-only job17700971 runs after any numerical exit. Remote batch
is the path below plus `_terminal_10`; local launch contracts are in
`terminal_10/`. New sources are `code/cluster/run_e5f_fixed_terminal_horizon.py`
and `code/model/tools/build_e5f_fixed_terminal_horizon_figures.py`.
No result from this additional arm is yet accepted. The original100-period
transition continues. Supplemental plots require visual QA after collection;
Python compilation passed before launch.

**18:22 UTC update:** replacement17699174 has verified the new stationary
endpoint, including all one-step population, queue, household, market and fiscal
checks. Its first saved-point replay is exact; stationary-root continuation
took210seconds. `terminal_restart_v1/verified_endpoint_receipt.json` is the
authoritative result. The anchored100-period transition is now running after
its first valid full mapping. All four other jobs remain running. The short
historical fit has completed one of four windows; no shocked path is yet
accepted. This endpoint result establishes neither uniqueness nor convergence
of the full transition.

Author approved launch on September 13, 2026. Smoke job **17697809 passed** in
4m55s and automatically dispatched all five jobs. At 17:46 UTC, finite IRFs
17697888/17697889 and historical fits17697891/17697892 remain running. No shocked
path is yet accepted. Existing presentation results remain frozen.

The endpoint arm17697890 exhausted24 evaluations after10m39s while improving:
maximum unscaled residual1.126e-5 versus5e-8 required. Replacement **17699174**
resumes the saved best coordinates and Jacobian, requiring an exact native
first-mapping replay before continuing. It retains the original one-hour
endpoint deadline18:12:42UTC and six-hour arm deadline23:12:42UTC. Source:
`code/cluster/resume_e5f_original_queue_terminal.py`. Local receipts are under
`terminal_restart_v1/`; remote replacement is the batch path below plus
`_terminal_restart_v1`. The original pinned solver and other four jobs are
unchanged. The earlier `prior_root_receipt.json` remains failed; only the new
`verified_endpoint_receipt.json` certifies the endpoint. Do not redispatch the
arm while its replacement transition is active.

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

Python compilation and native smoke passed. `smoke/summary.json` records
relative no-shock distribution drift4.84e-8 and adjusted queue drift4.41e-7.
The terminal solver reproduces the original stationary coordinates within
1.59e-6 relative. The shocked endpoint is verified; shocked transitions remain pending.
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
