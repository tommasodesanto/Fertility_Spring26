# Isolated current Estate-A transition on Torch

Current deployment defaults are v6 and use the twelve-start panel below. The
v5 launcher is frozen and its native smoke remains the reference for the exact
unchanged native identity and plans. Earlier standalone steps below document
that reference workflow; v6 uses `panel_launch_torch.sh`.

This launcher runs the current one-birth Estate-A adapter through the retained
one-permanent-shock controller. It creates one job at a time in a new immutable
root, `/scratch/td2248/projects/current_estate_transition_20261003_v5`. It does
not mutate or cancel existing calibration arrays or reuse historical transition
controller states or Jacobians.

The build uses the tested Estate-A calibration attempt3 archive with SHA-256
`5cb99a84f1aa51dab49462292ab31756f9183776e30d2d327844ae2af2c5fc44`.
It verifies every parent file, overlays the complete current birth-count and
transition-readiness packages, retains recovered historical dependencies,
includes the saved one-birth net-estate case and recursively includes declared
plan/handoff pins. Cache files are excluded. The builder materializes the exact archive payload locally and verifies it
before issuing the preparation receipt. The inventory pins every source,
data file, plan and launcher; host and mounted-container checks must both pass.

Preparation is deliberately separate from staging and submission. From the
repository root, with final, frozen smoke and fit plan paths:

```sh
python3 code/cluster/estate_birth_transition/deploy.py build \
  --smoke-plan /absolute/path/smoke_plan.json --fit-plan /absolute/path/fit_plan.json
python3 code/cluster/estate_birth_transition/deploy.py verify
python3 code/cluster/estate_birth_transition/deploy.py stage
ssh torch /scratch/td2248/projects/current_estate_transition_20261003_v5/launch_torch.sh preflight
python3 code/cluster/estate_birth_transition/deploy.py submit --mode smoke --wall-seconds 5400
```

`stage` writes the new remote directory but never submits a job. `preflight`
verifies the pinned smoke plan with zero model solves. Actual smoke invokes the
same current driver, fresh 12-date measured seed and six-date native loop with
six path maps; it does not start fit automatically. Inspect its actual receipt,
root/accounting/replay evidence and launcher terminal receipt before proceeding.
Copy `results/smoke/run/smoke_receipt.json` locally, then the lead may record an
explicit review and submit the separate bounded fit:

```sh
python3 code/cluster/estate_birth_transition/deploy.py review-smoke \
  --smoke-receipt /absolute/path/collected_smoke_receipt.json \
  --gate /absolute/path/fit_review_gate.json --reviewer lead
python3 code/cluster/estate_birth_transition/deploy.py submit \
  --mode fit --wall-seconds 21600 --gate /absolute/path/fit_review_gate.json
```

The review gate binds the exact smoke receipt, stage inventory, both plans and
current runtime identity. The remote fit launcher also compares the gate with
the actual remote smoke receipt. Fit cannot launch against a changed archive or
source identity. Finite-horizon diagnostic fit is not a certified production
policy result; the controller retains that distinction.

Each job uses eight CPUs, at most 96 GiB, eight Numba threads with single-threaded BLAS/OpenMP and a
strict external timeout. Submission requires an explicit wall budget of at most
six hours, including at least 120 seconds beyond the controller's total budget.
Plans must also pin individual stage budgets and iteration/policy-call caps.
Start/terminal receipts record job ID, mode and deadline; launcher heartbeat is
written every five minutes. The unchanged controller writes latest completed
and best-so-far summaries as evaluations complete. Existing result directories
are refused, and atomic local submission receipts prevent duplicate submission,
including automatic retry after an ambiguous or failed `sbatch` response.

No final bundle should be built while adapter, driver or either plan is still
being edited. A changed source requires a newly reviewed bundle and smoke,
not a relaxed source pin.

## Launcher-only revision from preserved v4

The failed v1 and v2 preflight roots is retained. For unchanged frozen sources and plans,
prepare a small v5 launcher delta locally instead of uploading the saved case
again:

```sh
python3 code/cluster/estate_birth_transition/deploy.py revise-launcher \
  --from-stage output/model/transition_readiness_v1/current_baseline_20261003/deployment_v4
python3 code/cluster/estate_birth_transition/deploy.py stage
```

Revision hardlinks immutable local source inputs, writes new launcher/deployer
and inventory bytes, and verifies all 697 source/data pins. The delta archive
contains only those three revised files. Staging first authenticates the v4
inventory and source payload, then copies its source tree on Torch with
`cp -a --reflink=auto` into the new absent v5 root and overlays the small delta.
No old run outputs or controller states are copied.

The consolidated bind manifest mounts explicitly complete pinned code packages, the current saved-case
subtree and current plan directory, with supplemental uncovered files afterward. It proves the effective
source origin of every original file pin and refuses duplicate targets or source
bind arguments above 100 KB (with 10 KB reserved for fixed binds). Historical output files remain individual overlays so frozen authenticated
siblings remain visible. The complete bind list stays below the 128 KiB
per-environment-string limit. The fixed publication inputs are mounted last, and any
staged pin within that directory must agree with the external input bytes.
Inherited Apptainer/Singularity bind environment variables are cleared before
execution. Host and container verification still hash every staged file; an
actual native import/preflight remains required before smoke submission.

Copied deployers resolve the canonical mounted repository root from adjacent
`inventory.json`; they do not assume `/work/deployment` has repository depth.

An optional `revise-launcher --fit-plan /absolute/path/fit_plan.json` replacement
changes only the pinned fit-evaluation and native-call caps. All other plan
fields, runtime identity, smoke plan and source/data bytes must match. The fit
file is detached from local hardlinks before writing; its file/plan SHA pins are
updated, and only that file is added to the launcher delta. This does not modify
the previous stage or the author-owned original plan.

The 96 GiB allocation requires Torch partition `cl`; `cs` rejected it before creating a job. Submission now explicitly selects `cl`, including when using the previously frozen v4 launcher. Smoke 19138166 was submitted with that explicit override. Read current receipts before any further submission.

## Frozen-source revision for the eight-thread native runtime

After the current runtime, handoff and both plans are frozen, `build --from-stage`
creates a source delta in the default `deployment_v5` location. Unlike the
launcher-only revision, it verifies and includes refreshed runtime/source and
handoff pins from both new plans:

```sh
python3 code/cluster/estate_birth_transition/deploy.py build \
  --from-stage output/model/transition_readiness_v1/current_baseline_20261003/deployment_v4 \
  --smoke-plan /absolute/path/frozen_smoke_plan.json \
  --fit-plan /absolute/path/frozen_fit_plan.json
python3 code/cluster/estate_birth_transition/deploy.py verify
python3 code/cluster/estate_birth_transition/deploy.py stage
ssh torch /scratch/td2248/projects/current_estate_transition_20261003_v5/launch_torch.sh preflight
```

The old payload is authenticated, unchanged local inputs use hardlinks, and
changed files are detached before replacement. The delta archive contains only
changed source/data/plan files and new launcher/deployer/inventory bytes. Remote
staging authenticates and copies the preserved v4 source tree on Torch, then
applies that delta. Every resulting file still passes host/container verification.
This avoids another upload of the approximately 229 MB saved-case payload.

The v5 job allocates exactly eight CPUs in partition `cl` and exports
`NUMBA_NUM_THREADS=8`; BLAS and OpenMP stay at one. Zero-solve preflight also
checks Numba's actual thread count. It remains one job running one controller,
with no array or simultaneous model solves. The existing deadline, stage/iteration
caps and all scientific gates are retained. An actual v5 smoke and explicit lead
review are required before v5 fit.

## v6: twelve independent scalar fits

v6 runs twelve scalar fits from distinct `fit_start_psi` values using the original
controller's `run()` method. The author-selected ratios to baseline
`psi_child=0.17892072066041628` are
`[0.1,0.3,0.5,0.65,0.75,0.82,0.88,0.94,1.0,1.1,1.3,1.6]`.
Each start has both strict diagnostic horizons 24/32, the same bounds, targets,
gates, 12-evaluation cap, 20,000-native-call cap and 21,480-second controller
budget inside a six-hour job. It independently reconstructs the reference and
measured seed; there is no shared mutable controller state.

The panel config schema is `estate_a_transition_panel_v1`. It contains `plan`
(an absolute fit-plan path and SHA), unchanged native `identity`, `panel_source`
(the absolute new panel-driver path and SHA), and exactly twelve ordered
`guesses` records with `index` and distinct `psi`. The native handoff and both
plans remain byte-exact v5 inputs. Panel orchestration is separately pinned in
the config and stage inventory.

After freezing the config and panel driver:

```sh
python3 code/cluster/estate_birth_transition/deploy.py build \
  --from-stage output/model/transition_readiness_v1/current_baseline_20261003/deployment_v5 \
  --smoke-plan /absolute/path/unchanged_smoke_plan.json \
  --fit-plan /absolute/path/unchanged_fit_plan.json \
  --panel-config /absolute/path/panel_config.json
python3 code/cluster/estate_birth_transition/deploy.py verify
python3 code/cluster/estate_birth_transition/deploy.py verify-panel
python3 code/cluster/estate_birth_transition/deploy.py stage
ssh torch /scratch/td2248/projects/current_estate_transition_20261003_v6/panel_launch_torch.sh preflight
```

Review the actual v5 smoke and panel mock tests before recording a v6 gate. The
panel test receipt must contain `status: "PASS"`, the exact
`panel_source_sha256`, `controller_run_reused: true`, and
`effective_plan_only_added_fit_start_psi: true`. The bridge requires identical
native identity, both plan pins and every existing source/data pin; any change
requires a new native smoke.

```sh
python3 code/cluster/estate_birth_transition/deploy.py review-smoke \
  --smoke-receipt /absolute/path/actual_v5_smoke_receipt.json \
  --native-stage output/model/transition_readiness_v1/current_baseline_20261003/deployment_v5 \
  --panel-test-receipt /absolute/path/panel_mock_test_receipt.json \
  --gate /absolute/path/v6_panel_review_gate.json --reviewer lead
python3 code/cluster/estate_birth_transition/deploy.py submit \
  --mode panel --wall-seconds 21600 --gate /absolute/path/v6_panel_review_gate.json
```

Submission uses `--array=0-11%12 --nodes=1 --ntasks=1`: at most twelve nodes,
each with eight CPUs and 96 GiB. This is below the user's 48-node cap. Each task
has its own `results/panel_guess_INDEX`, start/terminal receipts, five-minute
heartbeat and original controller checkpoints. The atomic panel submission
receipt refuses duplicates or automatic retries. The remote gate authenticates
the actual v5 smoke receipt and v5 inventory; the panel launcher additionally
verifies all v6 source/data, config and native-plan pins.

The driver collector can combine saved results with
`transition_panel.py --collect --inputs RECEIPT_OR_DIRECTORY ... --output OUT`.
It rejects mixed contracts and ranks accepted diagnostic fits by their final
window gap. No panel result is promoted automatically to a production policy.
