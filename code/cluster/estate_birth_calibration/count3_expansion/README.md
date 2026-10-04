# Additional count-three Estate-A searches

**Blocked October 4, 2026, about 00:48 New York.** The 412-file remote source
check and local zero-solve seed/collector tests passed, but all ten parent
tasks failed before final native verification. Controller **19139361** remains
PENDING on **19136605**, whose parent dependency cannot be satisfied. The new
count-three smoke has not run,
and there is no five-task production job ID or calibrated result. The
[stage](../../../../output/model/experiments/birth_count_choice/estate_a_count3_expansion_20261004_v1/deployment/stage_receipt.json),
[storage](../../../../output/model/experiments/birth_count_choice/estate_a_count3_expansion_20261004_v1/deployment/staging_storage.json),
[controller](../../../../output/model/experiments/birth_count_choice/estate_a_count3_expansion_20261004_v1/deployment/controller_submission.json),
and [validation](../../../../output/model/experiments/birth_count_choice/estate_a_count3_expansion_20261004_v1/deployment/validation_receipt.json)
receipts pin this pending stage.

This isolated stage adds five count-three optimizer chains to the five count-three
and five one-birth chains already planned in the overnight continuation. It
depends on successful controller **19136605**. That controller must first
authenticate all ten parent endpoints and pass fresh, exact two-call native
smokes in both arms. This stage verifies those receipts, constructs its five
starts, runs one more two-call count-three native smoke, and conditionally
submits a five-task array. The existing ten-task production array can run
concurrently: at most **15 Estate-A production cores**, plus this stage's
one-core controller before its five production tasks begin.

There is **no new economic change** relative to the Estate-A count-three arm:
the same post-interest timing, estate definition, preferences, entry and income
distributions, transfers, floors, population-one closure, empirical target
values and numerical weights apply. All ten calibration coordinates remain
free with the same bounds, including annual beta in `[0.93,0.99]`. The larger
household choice set does not mathematically nest the old aggregate fit; these
starts explore optimization, not identification or an adopted specification.

The five deterministic starts use the *lowest native-verified endpoint in each
parent arm*, with native loss compared only within the same new target contract:

1. Best binary endpoint, evaluated under the count-three arm.
2. Coordinatewise midpoint of best binary and best count-three endpoints.
3. Best count-three endpoint with beta set to `0.94` and first-birth fixed cost halved.
4. Best count-three endpoint with `kappa_fert` and `kappa_fert_continuation` tripled.
5. Best count-three endpoint with `theta0` doubled, `psi_child` multiplied by `1.5`, and `h_P` reduced by `0.25`.

Every value is clipped to the unchanged bound, and preparation fails if a
generated start duplicates another new start or one of the five original
count-three endpoints. These seed transforms do not fix parameters during the
search. Each chain starts a fresh Nelder-Mead simplex; no optimizer state is
claimed to resume. The plan pins source task IDs, parent receipt hashes, the
base controller's plan and smoke hashes, all ten free bounds, and complete
target and weight fingerprints.

The derived driver is byte-identical to the reviewed continuation driver
outside its new `checked_plan` gate. It retains the full native GE objective,
500-call limit, 32-lifecycle cap per GE, 1,800-second final native reserve,
fresh selected-point check, exact native repeat, 14-row target fit, 31-row
parameter record and 17 standard plots. Every task stops before **October 4,
2026, 10:00 New York** (epoch `1791122400`); a late start with under 45 minutes
remaining fails before solving. Production is released only after the own-arm
two-call native smoke passes. No automatic retry, fallback, pruning, parameter
promotion or scientific adoption is authorized.

The old ten chains plus these five have an upper planning limit of 7,500
retained one-GiB cases. The expansion controller requires at least 7,800 GiB
free at submission; individual evaluators retain the 350 GiB shared free-space
floor. At the parent array's observed 3–5 minutes per full GE, the 2,500-call
extra cap represents roughly 125–208 serial core-hours, so the clock is likely
to bind first. The stage records actual scratch free GiB when staged.

The reviewed `stage_torch.sh` created the isolated remote source stage and
passed a zero-solve 412-file host pin check. The reviewed
`submit_controller.sh` queued one controller with `afterok:19136605`.
The remote stage path is
`/scratch/td2248/projects/estate_birth_count3_expansion_20261004_v1`.
Read `staging_storage.json`, `control/controller_submission.json`,
`control/start_plan_receipt.json`, `control/smoke_gate.json`,
`control/production_submission.json`, and `control/controller_terminal.json`
as they become available. Five result folders provide heartbeats, latest and
best summaries, complete native reports or failures, and launcher terminal
receipts. `collect_torch.py --stage REMOTE --mode production --out NEW_JSON`
checks one arm's plan, target, source and fresh native receipts without
changing calibration status.

Local zero-solve checks: `/opt/anaconda3/bin/python code/cluster/estate_birth_calibration/count3_expansion/test_expansion.py`, `bash -n code/cluster/estate_birth_calibration/count3_expansion/*.sh`, and `/opt/anaconda3/bin/python -m py_compile code/cluster/estate_birth_calibration/count3_expansion/*.py`.
