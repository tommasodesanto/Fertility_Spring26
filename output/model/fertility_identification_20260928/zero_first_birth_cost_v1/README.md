# Zero first-birth fixed-cost diagnostic

## September 29 launch

Torch exact-loop tests **18817270** passed all five cases; preparation
**18817279** passed without model solves. Calibration **18817312** is launched,
with saved-output readout **18817376** dependent on its termination. The readout
skips if no normalized zero-cost center exists. Configuration SHA256 is
`7035e2c15641a6cf3dd5ddecd3eac1094bbaabe54057a7cf788c72081d73ce19`;
its hard end is epoch `1790726504.582713`. No result is claimed at launch.
Monitor `monitor-zero-first-birth-cost-test` checks once every 30 minutes,
remains quiet for healthy unchanged progress, and pauses after both jobs end.
It cannot repair, restart, extend or promote the experiment.

This isolated experiment tests the full target fit of the original one-birth
model when `first_birth_fixed_cost` is fixed at exactly zero. It is
an **experimental restriction**, not an adopted calibration or evidence of
global infeasibility. The comparison anchor is original selected
`one_birth_024_gn1_0`, loss 7.826226594410982, fixed cost
0.35270914196085973, normalized `psi_child` 0.12184551693359474. The source
case and checkpoint stay on Torch at
`/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project/output/model/fertility_identification_20260928/two_stream_overnight_v1/run_v1/one_birth/one_birth_024_gn1_0/case`.

## Economic contract and identification

The **only proposed economic restriction** relative to that selected candidate
is `first_birth_fixed_cost=0`. The original one-birth opportunity, earnings,
initial wealth and income distributions, timing, transfers/floors, preferences
other than this cost, housing and credit primitives, and the original ten scored
moment definitions and weights remain in place. The original three validation
rows remain validation only. `psi_child` is separately normalized to completed
fertility 2.1 at **every** case, and the demographic renewal gate is binding.
Nine other coordinates may subsequently move within their original bounds;
those movements are recalibration, not a second economic specification change.
There are ten scored targets for nine free parameters. Numerical Jacobian rank
will be reported, but neither numerical rank nor target count establishes
statistical identification. The **2007 stationary reference — block0506, September 28 verified export** and dated
transition are untouched.

The native evaluator still receives its full ten-coordinate point and full
original parameter-bounds contract. The restricted controller carries nine
search coordinates and adds cost zero before the native binding and evaluation.
It checks the bound parameter object, native receipt, 31-row saved parameter
table and scientific identity for the requested cost. The native runtime saves
the bound parameter object in the checkpoint packet and builds
`parameters.csv` from that same object (`e5f_evening_calibration_runtime.py`,
native evaluation/export path). This checks saved model-state provenance without
downloading or separately hashing a large checkpoint. The anchor replay is
explicitly marked as the unchanged positive-cost control.

## Exact staged loop and budgets

One serial Torch worker uses 1 CPU, 24 GB, a 4-hour internal budget and a
4h10m Slurm wall allocation. Each case has at most eight stationary solves
and 35 minutes. There are at most **16** cases: one authenticated selected
anchor replay, one zero-cost evaluation holding the other nine coordinates
fixed, nine fresh finite-difference columns at that zero-cost center, three
damped Gauss–Newton proposals, and two final repeats. The original positive
cost control does not compete for selection. A zero-cost center must pass the
2.1 normalization and renewal gates before any derivative or refit starts.
The anchor replay must match the saved selected loss within 0.05 and each
weighted residual within 0.01. The same repeat thresholds apply at the end.

At the prior median of about nine minutes per case, sixteen cases take roughly
2h24m plus setup and export; the hard case ceiling would exceed the overall
budget, so the four-hour clock takes precedence. The final 70 minutes are
reserved for two bounded repeats. If time, resource caps or scientific gates
interrupt a stage, the controller reports an incomplete packet and does not
impute derivatives, restart, relax gates, or claim infeasibility. The
zero-cost center can be censored under the unchanged eight-solve/35-minute
cap; that outcome requires a separate decision before any restricted refit.

Every successful case keeps all 14 target-fit rows, 31 parameter rows and the
native 17 standard plots. `heartbeat.json` updates every 15 seconds while a
case runs. `latest_completed.json` and `best_so_far.json` update after each
case; `first_stage.json`, `jacobian_0.json`, `selected_before_repeats.json`
and `FINAL.json` appear as their stages complete. The selected checkpoint hash
is checked once on Torch before repeats; large checkpoint files stay on Torch.
No restart or automatic extension is implemented.

## Torch preparation and launch sequence

Run these **after lead review** from the Torch stage root, using the repository
`code/cluster/torch.sh` workflow for sync and submission. `run.sh` refuses a
local invocation. Tests use synthetic subprocesses and no model solves. Then
prepare the exact configuration and inspect its receipt before submitting the
single search job. Preparation is exclusive-create and pins the prior search
configuration, original evaluator/source manifest, full target contract,
selected small artifacts and new controller/worker files. Neither preparation
nor tests imports the model.

```sh
sbatch --output=zero_cost_tests.log output/model/fertility_identification_20260928/zero_first_birth_cost_v1/run.sh tests
sbatch --output=zero_cost_prepare.log output/model/fertility_identification_20260928/zero_first_birth_cost_v1/run.sh prepare
# After both PASS and review of preparation_receipt.json:
EXPECTED_ZERO_COST_CONFIG_SHA256=<exact preparation_receipt config_sha256> \
  sbatch --export=ALL,EXPECTED_ZERO_COST_CONFIG_SHA256 \
  --output=zero_cost_search.log \
  output/model/fertility_identification_20260928/zero_first_birth_cost_v1/run.sh search
```

The preparation hard end is four hours after preparation. Submit promptly;
otherwise create a newly reviewed preparation rather than editing an expired
configuration. These commands describe the reviewed launch procedure; do not
repeat them to restart the already launched job.

## Saved-output readout

After `zero_cost_center` succeeds, run `postprocess.py` **directly from the
Torch `/scratch` path outside Apptainer** with one CPU, 4 GB and a ten-minute
limit. Re-run it after `FINAL.json` exists. It imports only the saved-data
projection code from the earlier comparison, verifies that source file's
pinned SHA, and reads authenticated small saved tables and observers.
It never opens a checkpoint, imports the model, or changes the native 17 plots.
The output under `readout_v1/` contains the complete 14-row fit and 31-row
parameter comparison, the four-panel fertility figure, other lifecycle
profiles, plotted CSVs and `verification.json`. Until two final repeats pass,
the selected result is labeled **provisional**. A failed zero-cost center has
no readout and must be described from `first_stage.json`/case failure evidence,
without an infeasibility claim.

```sh
sbatch --account=torch_pr_570_general --partition=cs --cpus-per-task=1 \
  --mem=4G --time=00:10:00 \
  --wrap='module load anaconda3/2025.06; /share/apps/anaconda3/2025.06/bin/python /scratch/td2248/projects/fertility_night_calibration_20260928_v1/project/output/model/fertility_identification_20260928/zero_first_birth_cost_v1/postprocess.py'
```
