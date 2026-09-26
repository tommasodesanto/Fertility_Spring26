# Four-arm utility comparison

The author delegated tonight's experimental choice; the reviewed fixed-rent
utility normalization is chosen for the comparison. The design and assumptions are in
[`docs/model/utility_four_arm_preparation.md`](../../../docs/model/utility_four_arm_preparation.md).

## Active overnight run

Torch array **18567879** was submitted once. All four arms entered their first
exact smoke on September 26 at 00:13 EDT. The shared search cutoff is 06:43,
the repetition cutoff 07:58, and the absolute finish 08:13 EDT. The source is
commit `d5dbf04d68ff000e1a1b8e66994cdffd31a01580` on
`codex/utility-four-arm-preparation`. The immutable approved contract hash is
`c3fa4d5b7925a54a747c72511538b8182029e1493e4d02e5ef703a30fcecc28a`.
`launch_v1/` retains the approved contract, freeze receipt, submission receipt,
and shared clock. The four slot assignments are `floor_linear`,
`floor_concave`, `shares_linear`, and `shares_concave`, respectively.

At this recorded launch observation, all four controllers had fresh heartbeats
and one active smoke; no objective had yet completed. This is launch evidence,
not a successful calibration or scientific-repeat certification. The common
smoke barrier gates the full search. Each arm reserves ten CPUs and 120 GiB;
the cluster independently enforces deadlines and writes latest/best summaries,
both selected repetitions and reports. No automatic resubmission is allowed.

This task owns the half-hourly follow-up
`check-overnight-utility-comparison`. It stays quiet on unchanged state, reports
material completion/failure, and stops after the final readout or a bounded
failure report by 09:00 EDT. The main research task owns canonical status
integration. Actual-result visual review is still pending.

## Preparation and reporting evidence

`preparation/` retains the first adapter-only check. The integrated
`preparation_v2/validation_receipt.json` records Torch job `18567413`: 54 tests,
four native arm preflights, two structural seeds per arm, zero equilibrium
solves. `preparation_v2/budget.json` records the eight-hour, 652-objective
finite design. Its `contract.json` is the checked preparation snapshot and
deliberately refuses objectives; a separate approved launch contract is
required. Subsequent source edits require a new snapshot, preserving old pins.

Remote preparation root:
`/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/utility_four_arm_preparation_20260925_v2/`.
The current checked packet is `checks_v2/`. Large checkpoints remain remote.
The approved launch is `launch_v1/`; active results are `results/run_001/`
beneath this same remote root.

`pdf_layout_validation/` contains a visibly labeled historical layout fixture,
its rendered pages and dependency receipt. This packet is not a new fit or
calibration result. The first fixture attempt (`18567452`) stopped because
ReportLab was missing; its partial outputs are retained. Only the reporting
environment is corrected before a fresh layout check.

The fresh reporting-only job `18567721` passed in 15 seconds, using ReportLab
5.0.1. Its final `visual_review_receipt.json` records all 22 pages, 13 target
rows, 29 parameter rows and 17 byte-identical historical figure images. New
tables and page layout have no clipping or overlap. Crowded income-state
legends and overlapping income-state labels are inherited plot limitations;
the stable figure set was preserved. This validates report pagination, not
the new economic results.

All targets, weights and measurement definitions inherit the immutable
September25 nightpair. The common adopted fiscal change and experimental
preference changes are separately disclosed in the design. Pending estate,
entry-wealth and mortality proposals are not inserted into these arms.

## Regenerate a completed arm's full diagnostic packet

Run on Torch after a selected original and both repeats exist. The collector
rechecks the saved scientific objects and all tables, then renders the same
17 figures and every PDF page without a new equilibrium solve. Use a fresh
output directory; it refuses to overwrite an existing report. For example:

```bash
module load anaconda3/2025.06
utility_root=/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/utility_four_arm_preparation_20260925_v2
utility_reference=/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/nightpair_20260925_v1
export PYTHONPATH="$utility_root/tools:$utility_reference/source/code/model/tools:$utility_reference/source/code/model:/scratch/td2248/commute_pdf_qa_deps"
export EXPECTED_UTILITY_COMPARISON_CONTRACT_SHA256=c3fa4d5b7925a54a747c72511538b8182029e1493e4d02e5ef703a30fcecc28a
python3 "$utility_root/tools/collect_e5f_utility_comparison.py" --stage collect --contract "$utility_root/launch_v1/contract.json" --arm floor_linear --run-root "$utility_root/results/run_001" --output "$utility_root/results/run_001/floor_linear/review_regenerated_001"
```
