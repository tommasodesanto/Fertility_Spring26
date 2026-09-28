# Transition preparation for the frozen September 28 reference

Owner: this chat. The economic-analysis chat owns saved-state interpretation,
fixed-price housing-cost experiments and the +10% supply comparison. The
calibration chat owns calibration improvement. No shared active model file is
edited by this preparation. The September 14 tag/reference is read-only.

**Reference:** 2007 stationary reference — block0506, September 28 verified export.
The authoritative identity is
[`../fertility_identification_20260928/fixed_reference_manifest.json`](../fertility_identification_20260928/fixed_reference_manifest.json),
SHA256 `147f9e2cb20f66350f1ceaa16cb41f822041ec869676ef5d5b9d04f16e4190d4`.
The selected checkpoint is `repeat_0212_primary`, SHA256
`b15ba92dc60e3d5590d2beb6e05d36f71d17b20b1a432edc2c2db926a217309d`.
Keep this reference even if another calibration improves.

## Current deliverable

[`preparation_v1/transition_readiness.md`](preparation_v1/transition_readiness.md)
is the substantive readiness note. Its initial numerical check is deliberately
limited to reference replay and no-shock operators. It does not authorize an
old credit mode, transport old shocks or certify an equilibrium transition.
The final [six-page PDF](../../pdf/fixed_reference_transition_readiness.pdf)
contains the complete reference fit and parameter tables. Exact-control reuse
and the one-date no-shock map pass. The longer demographic and fixed-stock
maps, natural-credit implementation, and a policy endpoint remain unverified.

Author decisions on September 28: hold physical housing stock fixed while
prices/rents clear; remove artificial borrowing/down-payment limits while
retaining lifetime solvency and repayment. The historical shock treatment is
separate and pending. All preferences, including saved child benefit, remain
fixed in credit counterfactuals. No recalibration or measurement change.

## Reproduction and evidence

The isolated Torch overlay is
`/scratch/td2248/projects/fixed_reference_transition_20260928/preparation_v1`.
Its launcher binds the existing frozen project read-only and this small overlay
read-write at its matching output path. Nothing is copied from or written over
the reference checkpoint. Only compact receipts and the report return locally.

- `preparation_v1/run_preparation_v4.py` and `run_v4.sh`: final instrumented
  one-date, saved-credit no-shock map. One CPU, 24 GiB, six-minute process cap,
  two dated Bellman calls. No recalibration. Torch job **18739319 passes**.
- `preparation_v1/runs/18739319/suite_result.json` and
  `no_shock_1/{identity,dated_audits,rows,calendar_display_rows,policy_progress}.json`:
  complete compact receipts. The native machinery has a legacy 2023 display
  label; `calendar_display_rows.json` labels this diagnostic from 2007 without
  changing the numerical calculation.
- `preparation_v1/runs/18738593/reuse_control/result.json`: imported exact
  control from the economics chat, 113 exact arrays, all 14 target rows,
  31 parameter rows and 17 plot hashes. No additional stationary solve.
- Preserved prior attempts: 18737340 stopped on a serializer representation
  mismatch before solving; zero-solve 18737715 identified the problem.
  Its 764 MB diagnostic JSON remains remote and must not be downloaded.
  Version 2 / 18738048 authenticated all 260 fields and completed one stationary
  solve, then an inherited diagnostic path hit the read-only mount.
- Version 3 / 18738593 passed control import, then its two-date map timed out
  at 180 seconds before a dated audit. Six-date and fixed-stock stages were
  not reached. This is not an economic or feasibility rejection.
- Version 4 changed the diagnostic design to one instrumented date with a
  six-minute cap. Its two policy calls took 53.753 and 51.025 seconds; the map
  took 111.860 seconds, 142.013 including setup. No automatic retry follows.
- `preparation_v1/build_report.py` and `render.sh` generate and rasterize the
  report on Torch from the frozen manifest and saved receipts, with no model
  imports. All six pages are visually reviewed before delivery.
- `preparation_v1/historical_user_evidence.jsonl`: selected September 12-13
  user messages from *Review quantitative model*, identified using chat tools
  and a narrow session-index query before searching its one candidate file.
- `preparation_v1/historical_source_identity.json`: old source/empirical pins,
  extracted from the historical manifest. These identify old evidence; they
  do not assert a fresh hash replay of the old remote source tree.
- `preparation_v1/history_review.md` and `runtime_review.md`: bounded read-only
  worker reviews. The lead verified the relevant sources. Correction to the
  history review: an old **fixed housing-supply curve** still had elasticity
  0.630; it was not a fixed physical stock.

The full reference [14-row fit](../fertility_identification_20260928/resume_v1/selected_export/primary/target_fit.csv),
[31-row parameter table](../fertility_identification_20260928/resume_v1/selected_export/primary/parameters.csv)
and [17 standard plots](../fertility_identification_20260928/resume_v1/selected_export/primary/standard_diagnostics/)
remain unchanged. No claims of calibration optimality or policy effects follow
from this preparation.

## Observed checks and remaining gates

The one-date map has housing residual 1.760e-9, PAYGO residual 3.172e-13,
distribution L1 change 4.010e-14, and zero mass-accounting error and inherited
projection. All dated household, purchase/repayment, occupied-value,
probability and actual-next-cohort estate-funding gates pass. Its adjusted
entry queue changes by 2.445e-8, consistent with the retained tiny baseline
renewal discrepancy. These receipts validate a constant-path operator only.

Do not launch a production counterfactual from these scripts. Required next
steps are a six-date loop smoke spanning both entry lags; the fixed-stock
closure check; an explicit Boolean natural-solvency implementation and debt-grid
verification; a demographic/fiscal/housing endpoint at fixed preferences; then
finite-path clearing and horizon convergence. Historical fixed-stock work also
needs an agreed shock contract. The old long-horizon refit accepted zero shocks.
The provisional estate counterparty/physical settlement remains outstanding.

A pure supply-intercept increase has a scale-invariant endpoint candidate under
the present single-market closure. See the note for the checked linearity
argument, scale-breaking objects and remaining native scaling test. It is not
a transition, uniqueness or stability result. The economics chat owns that arm.

## Commands

One-date verification only (creates a new run folder keyed by Slurm job ID):

```sh
ssh torch 'sbatch --output=/scratch/td2248/projects/fixed_reference_transition_20260928/preparation_v1/slurm_%j.log /scratch/td2248/projects/fixed_reference_transition_20260928/preparation_v1/run_v4.sh'
```

Report regeneration from the settled receipts, without a solve:

```sh
ssh torch 'sbatch --output=/scratch/td2248/projects/fixed_reference_transition_20260928/preparation_v1/render_%j.log /scratch/td2248/projects/fixed_reference_transition_20260928/preparation_v1/render.sh'
```

All numerical work and rendering ran on Torch. No active model, empirical,
calibration, manuscript or September 14 reference source was edited. The saved
reference project stayed read-only. Only this packet and its final PDF are included in the bounded source-control
commit; automatic Git maintenance is disabled. Other existing changes in the
shared checkout are excluded.
