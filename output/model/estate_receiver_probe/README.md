# Estate recipient diagnostic

## Status — September 25, 2026

Author requested the three-way test discussed in the bequest-accounting task.
Torch access is restored. Native source/checkpoint preflight passed in job
18563842 (`results/v2/preflight`), and all nine focused fixture tests passed in
18564037 (`results/v2/smoke`). Comparison job 18564157 uses `results/v2/run`,
with exact control replay required before treatments. The control matched all
saved solution arrays and targets exactly; the valuation case also completed.
The recipient fixed point is running, so no full comparison is claimed yet.
The first preflight (18563806) stopped before
model import because Slurm relocated its script; the launcher now preserves the
submission project root explicitly. Model work, imports, tests, source hashing
and rendering stay on Torch. Existing reference results are untouched.

## Reference and economic changes

Reference: September 25 paired PAYGO experiment, `oasi_087510`, actual selected
case `worker_05/point_04`. Remote parent is
`/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/nightpair_20260925_v1`.
The control retains that experiment's tax rate 0.08751017424959717, all selected
structural parameters, fertility preference, B15 earnings, entrant joint wealth
distribution, wealth grid, housing supply, target rows/weights and numerical
gates. It does not implement the newly adopted pension ratio/tax rule and is
not an adopted calibration.

| Case | Bequest utility estate | Receipt by living households |
|---|---|---|
| control | `max(b' + q*h, 0)` | none |
| net_valuation | `max(b' + (1-psi)*q*h, 0)` | none |
| net_valuation_transfer | same net estate | equal anticipated period transfer to age cells starting in [45,65] |

Both economic changes are **experimental, not adopted**. The existing four-year
grid pays age-start cells 46, 50, 54, 58 and 62; this is not a precisely
age-integrated recipient profile. No family matching or random inheritance draw
is added. The existing current-income purchase rule includes the receipt in
`Y/R`, with no second credit to transaction wealth or next-period assets.

The receiving case solves `T*M = D`, where `M` is recipient household mass and
`D` is the period death-weighted positive estate flow net of selling costs.
Estates use post-saving wealth and the native post-transaction current measure;
terminal death has probability one. Entry wealth is not funded or changed by
this pool. There is no transfer at age 18.

The empirical bequest target remains the original **gross** observer in all
three cases. Net flows and receipts are supplemental accounting. Receipts do
not enter measured labor earnings. Thus a target-definition change cannot
masquerade as an economic response. The existing all-positive-estates versus
child-directed empirical-target mismatch remains disclosed and unchanged.

The fixed-preference comparison keeps normalized entrant/age composition and
reports the birth-based adult-entry replacement residual. Treatment cases are
not claimed to be closed-population stationary equilibria. It does not impose
fertility 2.1 by changing preferences. Policy transitions and any subsequent
recalibration remain a separate stage after this diagnostic.

## Implementation and verification

- `code/model/tools/e5f_estate_receiver_adapter.py`: isolated in-memory income
  and Bellman hooks, supplemental death accounting and bounded transfer loop.
- `code/model/tools/run_e5f_estate_receiver_probe.py`: pins the parent lock,
  source, targets and selected checkpoint; reuses its verified purchase runtime,
  observers, fiscal/market/household/purchase/mass/value gates and 17 figures.
- `code/cluster/submit_e5f_estate_receiver_probe.sh`: separate Torch stages.
- Two `test_*estate_receiver*.py` files test period/annual accounting, death
  timing, income exclusion, funded receipts, control nesting, loop/cap failures.
- `code/model/tools/build_e5f_estate_receiver_report.py` assembles a completed
  three-case packet into a PDF with all tables and 51 retained figures. It is
  source-reviewed; runtime rendering and visual QA await the completed run.

Source review corrected the worker draft's ancestor-versus-selected preference
comparison, remote paths, solve counter, unconditional transfer balance gate,
and confusion of mass conservation with demographic replacement. Source
preflight and all nine focused tests have now passed on Torch. The native control must reproduce
the selected solution arrays and all target model values exactly; failure stops
the experiment and does not authorize wider tolerances.

Budget: 20-minute preflight allocation, 20-minute smoke allocation (tests capped
at 10 minutes), then a three-hour run deadline with 42 native solves maximum:
one control, one valuation case, at most 40 receiver iterations. The reference
chosen equilibrium took 99.811 seconds, implying about 70 minutes at that speed
for 42 solves; nonlinear convergence, audits and graph export add uncertainty.
The earlier different-model receiver test took 54 minutes. No search or repeated
full-case battery is planned. Every iteration writes a receipt; a 55-second
heartbeat, latest completed case and descriptive best-so-far remain readable.
Nonconvergence, baseline mismatch or a numerical gate stops the run.

Completed cases will retain checkpoints, full 13-row target-fit and parameter
tables, estate accounts and the stable 17-figure packet. A readable comparison
PDF still needs assembly/visual QA after results exist; no numerical deliverable
is claimed from source preparation alone.

## Resume on Torch

Place only the four new Python files and
launcher in a **separate** task tree, preserving their `code/model/tools/` and
`code/cluster/` paths. Suggested remote task root:
`/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/estate_receiver_probe_20260925_v1`.
Never copy them into or modify the frozen `nightpair_20260925_v1/source` tree.

The launched comparison uses `results/v2/`; those directories must not be
reused. For a separately justified future run, choose a fresh version below.
From the separate task root, the stage sequence is:

```bash
export E5F_ESTATE_PROBE_STAGE=preflight
export E5F_ESTATE_PROBE_OUTPUT="$PWD/results/NEW_VERSION/preflight"
bash code/cluster/submit_e5f_estate_receiver_probe.sh --submit
```

Inspect completion, then:

```bash
export E5F_ESTATE_PROBE_PREFLIGHT_RECEIPT="$PWD/results/NEW_VERSION/preflight/preflight.json"
export E5F_ESTATE_PROBE_STAGE=smoke
export E5F_ESTATE_PROBE_OUTPUT="$PWD/results/NEW_VERSION/smoke"
bash code/cluster/submit_e5f_estate_receiver_probe.sh --submit
```

The smoke runs the exact driver loop and actual adapter fixed point with
deterministic fixtures; no full native equilibrium is solved in that stage.
After its success, the first production case is the exact native baseline
replay, and treatments may proceed only after it passes:

```bash
export E5F_ESTATE_PROBE_SMOKE_RECEIPT="$PWD/results/NEW_VERSION/smoke/complete.json"
export E5F_ESTATE_PROBE_STAGE=run
export E5F_ESTATE_PROBE_OUTPUT="$PWD/results/NEW_VERSION/run"
bash code/cluster/submit_e5f_estate_receiver_probe.sh --submit
```

Receipts pin the complete parent source/target/checkpoint and the new diagnostic
files. A source change invalidates smoke/preflight. Do not retry or enlarge the
budget without a changed, stated diagnosis. The launcher does not poll or create
a monitor; active job status is recorded above.
