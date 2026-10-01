The launcher reuses the calibrated frozen-root, grid-resolution base overlays, continuation source overlays, and authenticated input bundle. It adds only an inventoried transition source and verified calibration handoff overlay. It never writes the calibration deployment.

`stage.py` snapshots the transition source, handoff, verified selected packet, source manifest, and any prepared plan and annual/blocks inputs into `stage/`. Upload that directory to `/scratch/td2248/projects/transition_readiness_v1/current_floor/` with `rsync -az`. The inventory authenticates every uploaded source before container entry. Rebuild the snapshot after any source or plan change; do not mutate an inventory used by a running job.

Zero-model-call target-architecture preflight (run on the Torch login node):

```sh
bash /scratch/td2248/projects/transition_readiness_v1/current_floor/floor_launch.sh --mode preflight --seconds 300 --label native_import_new_label
```

After copying its `preflight.json` here as `native_import_preflight.json`, `make_plan.py` constructs the lead-specified diagnostic smoke plan. It takes the identity from the actual zero-call constructor, actual selected preference from the parameter table, and the retained complete annual target contract. It does not submit a job.

The lead must review and submit the final inventoried smoke plan. Pass `--time=01:31:00` and an explicit log output to `sbatch`, then use launcher arguments `--mode smoke --seconds 5400 --label smoke_v1 --plan /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/transition_readiness_v1/current_floor_preparation/deployment/smoke_plan.json`. Allocation is 4 CPUs/96 GiB; native numerical libraries remain one thread. The 90-minute actual-start timeout, per-stage controller limits, iteration limits, and 260 actual policy-call cap stop work without automatic retry. The launcher writes start/terminal receipts and a heartbeat every five minutes; native controller checkpoints preserve completed reference/seed inputs for explicit later reuse.

The short six-date smoke verifies the changed-shock equilibrium and fresh replay. Its terminal comparison remains diagnostic; it cannot certify the full 104/128-date transition or production readiness. No model job is authorized by running these preparation scripts alone.

Actual target-architecture import `native_import_v2` passed with zero policy calls and exit 0. Its receipts and log are retained in `native_import_v2/`; the verified constructor still requires native selected-state reconstruction. The generated `smoke_plan.json` passed the complete zero-call controller preflight. The final 89-file remote inventory was checked after upload. No Slurm job was submitted by this preparation worker.

Native checkpoints and measured-seed receipts use stable `/work/transition_runs/<label>` paths. All prior result folders are mounted read-only at `/work/transition_runs`, with only the current label mounted writable. `/work/results` remains a numerical cache alias. This lets a later explicitly pinned fit reuse completed readiness inputs without rewriting receipt paths.

After the first actual smoke stopped on a native utility export mismatch, the scoped runtime callback repair was deployed with updated controller failure counters. Actual `native_api_v1` passed all 23 required helpers with zero policy calls and exit 0 in 14.7 seconds. Its receipt and log are retained here. The regenerated 89-file inventory and smoke plan pin runtime `e1ecd57d…` and controller `6eafa7da…`; use a fresh `smoke_v2` label for the next lead-submitted model run. The six-part native identity and all numerical budgets remain unchanged.

Actual one-shock fit `actualfit_v1` was submitted as Torch job `18971392` after the exact mounted preflight authenticated the original native checkpoint, measured matrix, generation/consumer audit, and critical AST comparison. It uses horizons `[12,16]`, a five-hour actual-start limit, and 2,200 actual native calls. It began at 15:24:47 New York on October 1; the actual deadline is 20:24:47. Reference/seed restoration passed with zero new solves. The full 98-file executed snapshot is retained under `current_floor_runs/actualfit_v1/executed`. Candidate completion and a fitted shock remain to be observed.

Reproduce the compact read-only collector without importing the model or copying packet arrays:

```sh
ssh -4 -oConnectTimeout=20 -oBatchMode=yes torch '/share/apps/anaconda3/2025.06/bin/python -' < /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/transition_readiness_v1/current_floor_preparation/deployment/actualfit_progress.py
```

It reports launcher start/termination, failure, current phase, completed candidates, latest/best measurements, and exact 2023 checkpoint receipt paths while recursively suppressing source identity maps and arrays. The last collected snapshot is `current_floor_runs/actualfit_v1/progress_snapshot.json`.

V1 was cancelled after its six-map H12 derivative root stopped at score 1.762519 despite shrinking residuals. The completed baseline and all compact mapping/endpoint receipts remain archived. The true completed-call lower bound is 101; the interrupted twelve-date map gives a conservative upper bound of 125. V2 charges 125 against the original 2,200-call budget and retains the original October 1 20:24:47 New York hard deadline. Its proposed diagnostic controls are twelve path evaluations, full path damping, and twelve scalar-fit evaluations; endpoint damping, economics, targets, and gates remain pinned. Submission requires the corrected approved compatibility receipt, exact mounted zero-call preflight, and final lead review. The launcher now accepts `--deadline-epoch` and computes its external timeout at actual start.

Deployment recovery: an erroneous mirror upload deleted `current_floor/results` on Torch before interruption. No numerical job was active. The exact tree was restored from the authenticated 15:46 EDT Torch filesystem snapshot, excluding replaceable caches, then supplemented with later compact v1 evidence from the local archive. `recovery_receipt.json` records 1,169 files / 3,296,771,402 bytes and matching native reference, measured Jacobian, original preparation, and baseline 2023 checkpoint hashes. The local baseline checkpoint is also independently retained. Future uploads use `upload_stage.sh`, without deletion and with explicit results protection.

V2 exact mounted preflight passed with zero native calls. `actualfit_v2_delta.json` records the complete plan changes and hard deadline; `actualfit_preflight_v2/` contains the actual terminal, log, and validation receipt. The complete 101-file intended source snapshot is retained under `current_floor_runs/actualfit_v2/executed`. `submit_actualfit_v2.sh` is prepared for lead-approved submission and computes remaining time against the original deadline; it has not been executed.

After final lead review, `actualfit_v2` was submitted exactly once as Torch job `18974228`. It actually started at 16:22:00 EDT, with the unchanged original 20:24:47 deadline and 2,075 remaining-call cap. Native reference/seed restoration passed with zero new solves and candidate 1 began at the original baseline preference. `actualfit_v2_first_phase.json` records this actual evidence; `actualfit_v2_progress.py` is the compact reproducible collector. No further source staging is permitted while this generation runs.

Never use `--delete` for this task directory. Use `upload_stage.sh`, which explicitly excludes generated results and Slurm logs. Preserve each executed source inventory before staging another version.
