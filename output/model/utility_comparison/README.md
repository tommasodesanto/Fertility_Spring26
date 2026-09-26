# Four-arm utility comparison

## Collected comparison

The [complete four-arm readout](run_001/comparison.md) now links every target,
parameter, actual report and qualified visual review. Normal cluster access
returned on September 26 and both share exports were retrieved. All four
selected originals match each of their two repeats exactly. Every report has
17 standard figures and 22 reviewed pages; floor parameter tables have 28 rows
and share tables have 29. The original [bounded failure readout](run_001/readout.md)
and its receipt remain the historical account of the earlier access block.

**The intended recalibration was not achieved:** all searches stopped in their
initial populations, no DE generation ran, and coverage differs across arms.
The run also retains the old ACS mean-room target, rather than the author's
prior AHS choice. The readout records the next-contract correction and an
unresolved upper-wealth housing-policy feature in the share packets. Exact
repetition does not establish a preferred specification or grid robustness.
`run_001/collection_completion_receipt.json` records the final collection checks.
The heartbeat remains paused; all four original jobs are now accounted for.

The author delegated tonight's experimental choice; the reviewed fixed-rent
utility normalization is chosen for the comparison. The design and assumptions are in
[`docs/model/utility_four_arm_preparation.md`](../../../docs/model/utility_four_arm_preparation.md).

The main research task subsequently requested a separate recovery design.
[`docs/model/utility_four_arm_recovery_v1.md`](../../../docs/model/utility_four_arm_recovery_v1.md)
contains the proposed policy, reuse conditions and a smaller **eight-hour
option** with at most forty new points per arm. `recovery_v1/eight_hour_option.json`
is the current budget proposal; the earlier 42-hour calculation is an
unapproved upper-bound illustration, superseded as the default. This is preparation only,
not an approved launch, retry or extension of `run_001`. Both proposals predate
the accepted-input reconciliation and cannot authorize reuse under different
targets or observers.

The bounded machine repair subsequently passed **60 focused tests on Torch**
and a zero-solve check of the actual frozen native exception class. The
[repair receipt](recovery_v1/repair_test_receipt.json) preserves the two earlier
failed test logs, final source fingerprints, review scope and remaining limits.
The opt-in controller continues after narrowly evidenced candidate rejection,
preserves owned timeouts as unscored, and verifies checkpoint resume without
replaying completed attempts. A batch containing timeouts still stops before
DE selection. A corrected scientific contract and a separate bounded native
check remain necessary before another long run.

## Outcome in progress: planned recalibration not achieved

By the 01:45 EDT check, all four searches had stopped during their initial
populations. **The planned full eight-hour recalibration was not achieved;
no differential-evolution generation ran.** The frozen policy stopped an
entire arm after its first unexpected feasibility failure or time limit.
Successful saved points are still useful for a conditional comparison, but
repeatability of those points does not repair unequal search coverage or
establish a preferred utility specification.

Both floor arms have completed their selected-point repetitions and reports.
Each original matches each of its two repetitions exactly, and each report
contains all 13 target rows, 28 parameters, 17 standard figures and 22 pages.
Their actual reports have passed visual review with inherited legend/tick
overlaps documented in `export/visual_review_receipt.json`. Read each PDF
with the early-stop and three-successful-point limitation below. The two
share arms were finishing already-running cases under the original limits;
their final counts, repetitions and exports were still pending at this check.
The 02:14 EDT heartbeat now shows both share arms running their two selected
repetitions. Their initial trials are finished; export and actual-report QA
remain pending. The exact heartbeat and raw count snapshot is retained in
`run_001/status_20260926_0614.json`.

At 02:21 EDT, the next status read failed SSH authentication; the local Torch
control socket was absent. No new cluster state was obtained, so this is a
collection-access problem, not evidence of a controller failure. The main
task was notified, credentials/settings were left untouched, and the existing
monitor remains active. `run_001/collection_access.json` records the bounded
diagnosis. Final share-report verification waits for normal access to return.

## Launch and stop evidence

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

At the 01:15 EDT check, all four arms had passed both exact smokes and entered
the search. Both housing-floor arms then stopped new search dispatch after
the same initial trial (`initial_0005`) raised `InfeasibleThetaError`. At age
34, positive mass in states without a feasible allocation exceeded the
unchanged threshold: approximately \(1.194\times10^{-12}>10^{-12}\). The
census includes renter households with zero liquid wealth, one child at home,
and a negative budget slack of approximately 0.024 model units. The small
mass does not establish that this violation is harmless. This error was not
one of the three prespecified inadmissible-proposal classes, so the frozen
controller stopped rather than changing its classification or retrying.

Each floor arm attempted ten new initial points: two succeeded, seven were
rejected by the named housing-equilibrium gate, and one failed as above. Each
retained the verified smoke seed, leaving 29 initial slots and all 120 DE
trials unrun. Their selected results are therefore best among only three
distinct successfully scored points; they cannot support a fair comparison
of fully recalibrated specifications against more extensively searched arms.
The floor controllers moved to the two required repetitions and export. Both
share arms were still running at this check. Full failed-trial parameter
vectors and native feasibility censuses are preserved in each floor arm's
`run_001/<arm>/initial_0005_record.json` and `initial_0005_failure.json`.
These failures do not provide a completed failed-trial
equilibrium or a validated estate-accounting result.

Both share arms subsequently stopped new dispatch after `initial_0012`
reached its 3,100-second objective cap. The saved wall times were 3,100.932
and 3,100.582 seconds for linear and concave benefits, respectively. The
tracebacks explicitly show `TimeoutError: objective wall budget exhausted`
from the runner's alarm, subsequently wrapped by Numba as `SystemError`.
Consequently, the immutable controller records count these as `failed`, not
`timed_out`. The postmortem identifies the deadline root cause separately;
this is not evidence of an independent kernel defect. The corresponding
records, failure receipts, case plans and complete short tracebacks are
preserved under `run_001/<arm>/initial_0012_*`. The interrupted objectives
had completed ten and eleven stationary normalization evaluations,
respectively. All 39 new initial points had already been dispatched in each
share arm; remaining siblings finish before the existing repeat/export path.

This task owns the half-hourly follow-up
`check-overnight-utility-comparison`. It stays quiet on unchanged state, reports
material completion/failure, and stops after the final readout or a bounded
failure report by 09:00 EDT. The main research task owns canonical status
integration. Actual floor-result visual review is complete; share-result
review remains pending.

The requested supplemental birth-versus-wait diagnostic is unavailable from
these saved outputs. The frozen solver computes the separate action values
in temporary `Vfa` and `V2` arrays, then saves inclusive values and choice
probabilities. It does not retain the required action-value arrays and
eligibility masks. A bounded source-schema review, independently checked by
the lead, is recorded in `run_001/supplemental_birth_wait/availability.json`.
No probabilities were inverted, no extra solve was run, and no zero-shock
equilibrium claim is made. Exact repetitions also do not establish grid
convergence. The separate borrowing and wealth-grid reviews remain distinct
from this stationary utility experiment.

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
