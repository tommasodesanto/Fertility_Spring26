# Runbook: overnight psi>=0 twin arrays (E5F floor arm + repaired-E5 control)

Execution-only handoff. All implementation is committed and lead-reviewed;
do not edit model code, launchers, targets, or weights. Your job is to
submit, health-check, and report. If any gate below fails, stop and report —
do not restart, broaden, or fix silently.

Author decision being implemented: psi_child (child utility flow) is
restricted to be nonnegative in both overnight arms. The mechanism question
(floor versus tilt) is settled tomorrow by comparing the two certified
winners under this common restriction.

Cluster: torch, account `torch_pr_570_general`, from the repo checkout at
`~/Fertility_Spring26` (scratch conventions as usual). First action:

```bash
cd ~/Fertility_Spring26 && git pull --ff-only
```

Confirm `git log --oneline -1` includes the psi-floor commit (message
"Add author psi floor restriction..."). The runs are invalid without it.

## Arms

Both arms: J=17, Nb=120, search 10/1e-4, tight winner 40/2.5e-5 x2 exact
repeats (built into the chain), 8 chains x 225 minutes, seeds from the
repaired-E5 certified winner
(`output/model/eqscale_seq_e5_maturation_repair_recalibration_20260805/report/results.json`).
Note: that winner has psi_child = -0.558; under the psi floor the seed clips
to psi = 0 and the start-mix spreads chains from there. This clipping is
expected, not an error.

- **Arm A (E5F child-room floor):** 9 free parameters (tilts out, one
  rooms-per-child floor in, domain [0.10, 1.80] log), psi_child in [0, 3].
- **Arm B (repaired-E5 control):** the unchanged 10-parameter tilt
  specification, psi_child in [0, 3]. This is the clean comparator; without
  it the floor-versus-tilt comparison would be contaminated by the psi
  restriction itself.

## Phase 1 — Torch smokes (submit both, ~20 min)

```bash
cd ~/Fertility_Spring26/code/cluster
sbatch --export=ALL,E5_PSI_MIN=0.0,E5F_RUN_TAG=psinneg_smoke_20260806,E5F_SEED_RECORD=$HOME/Fertility_Spring26/output/model/eqscale_seq_e5_maturation_repair_recalibration_20260805/report/results.json submit_intergen_e5f_child_room_floor_smoke.sh
sbatch --export=ALL,E5_PSI_MIN=0.0,E5_REPAIR_RUN_TAG=psinneg_smoke_20260806,E5_REPAIR_SEED_RECORD=$HOME/Fertility_Spring26/output/model/eqscale_seq_e5_maturation_repair_recalibration_20260805/report/results.json submit_intergen_e5_maturation_repair_smoke.sh
```

(If the repaired-arm smoke launcher has a different filename, use the smoke
launcher that ran on 2026-08-05 for job 15398044-era smokes; same env vars.)

**Pass gates, checked in each smoke chain's `metadata.json` and
`cases.jsonl`:**

1. `metadata.json` records `"psi_min": 0.0` and, for Arm A, the E5F profile
   block (9-parameter domain including `hbar_child_rooms`).
2. Every completed smoke case is strict (`n_strict == n_cases_completed`),
   no error log entries.
3. No case's theta contains `psi_child < 0`.

If any gate fails: report the failing file verbatim and stop.

## Phase 2 — Production arrays (only after both smokes pass)

```bash
cd ~/Fertility_Spring26/code/cluster
A=$(sbatch --parsable --export=ALL,E5_PSI_MIN=0.0,E5F_RUN_TAG=psinneg_20260806,E5F_SEED_RECORD=$HOME/Fertility_Spring26/output/model/eqscale_seq_e5_maturation_repair_recalibration_20260805/report/results.json submit_intergen_e5f_child_room_floor.sh)
sbatch --dependency=afterany:$A --export=ALL,E5_PSI_MIN=0.0,E5F_RUN_TAG=psinneg_20260806 submit_intergen_e5f_child_room_floor_collector.sh
B=$(sbatch --parsable --export=ALL,E5_PSI_MIN=0.0,E5_REPAIR_RUN_TAG=psinneg_20260806,E5_REPAIR_SEED_RECORD=$HOME/Fertility_Spring26/output/model/eqscale_seq_e5_maturation_repair_recalibration_20260805/report/results.json submit_intergen_e5_maturation_repair.sh)
sbatch --dependency=afterany:$B --export=ALL,E5_PSI_MIN=0.0,E5_REPAIR_RUN_TAG=psinneg_20260806 submit_intergen_e5_maturation_repair_collector.sh
```

Record all four job IDs in your report.

## Phase 3 — Health checks (twice: ~30 min in, and before you finish)

- `squeue` shows 8+8 chains running (then collectors pending on
  dependencies).
- Each chain's outdir gains `cases.jsonl` lines and a `latest_completed_case.json`
  heartbeat; if any chain writes nothing for 30 minutes, flag it (do not
  resubmit).
- Spot-check one early case per arm: theta respects `psi_child >= 0`; Arm A
  thetas contain `hbar_child_rooms` and no `delta_alpha` keys.

## Phase 4 — Morning report (report only; adopt nothing)

For each arm: chain eligibility count, winner tight loss + residual +
exact-repeat evidence, full 12-row target table, all free parameters with
bounds and near-bound flags (watch: psi_child at the 0 bound; Arm A's
hbar_child_rooms near 0.10 or 1.80; theta0/theta1 lower bounds). Do not run
policy packets, do not update CALIBRATION_STATUS beyond appending the run
record, do not modify the maintained model. The lead audits both reports and
the author decides between devices.

## Explicitly out of scope tonight

Any change to targets, weights, evaluator settings, launcher internals, or
model code; any E6 arm; any policy counterfactual; any comparison quoted
across different psi domains (last night's 249.19 winner is psi-unrestricted
and is NOT the comparator for tonight's arms).
