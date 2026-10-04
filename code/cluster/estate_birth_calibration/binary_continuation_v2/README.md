# One-birth Estate-A continuation v2

**Live launch, October 4, 12:39 New York:** exact two-call native smoke
`19163407` passed. Production array `19164416` (`0-19%20`) was submitted
with a five-minute `--begin` delay so its job-ID receipt was durable before
any task could start. All 20 tasks started and wrote launch, start-contract
and first full-GE heartbeat receipts. Tasks 14 and 19 then failed their first
native GE on the unchanged $10^{-12}$ dead-node mass gate; 18 remain running
as of 12:44 New York. Neither failed task was restarted. The first
held submission `19164358` was canceled before any task started because
Torch returned `Unspecified error` on `scontrol release`; both submission
receipts are retained in the local deployment folder and remote `control/`.
Do not run `submit_torch.sh` again. Read the current Slurm queue, the
`control/production_submission.json` receipt, and each
`results/production_binary_chain_N/run/{heartbeat.json,latest_completed.json,best_so_far.json}`
to monitor. Investigate any running task without a new heartbeat for 30 minutes.
The native `completed.json` and collector, not a best-so-far checkpoint,
decide whether the verified loss is below 13.

This package prepares a new search stage for the one-birth Estate-A calibration. It derives from the authenticated October 4 recovery stage (array `19141024`, inventory SHA `974a14e243da6a2ad0572bb9825b47ab349828f9144cbbad69e74b40a4408b22`). The objective remains the same 14-row target and weight contract, ten free parameters and bounds, selected-source SHA, native full-GE evaluator, 32-lifecycle cap and fresh selected-point exact-repeat gate. Economic inputs, target measurement, objective and solver do not change; the staged driver adds only plan authentication and search control.

The 20 deterministic starts contain five verified recovery winners, five feasible saved binary checkpoints from the failed parent array `19127370` (clearly labeled provisional and authenticated against its inventory, contracts and case logs), and ten bounded perturbations around verified recovery chain 1. The builder authenticates every source receipt, keeps every coordinate inside the unchanged bounds, rejects duplicate or too-close points and records all source hashes. Search starts from these vectors with fresh Nelder-Mead simplices; optimizer state is not resumed.

Each production task has one CPU, 24 GiB, a 12-hour wall clock beginning at task start, at most 500 objective calls, an 1,800-second native reserve and no automatic retry or extension. Before release, the submitter requires at least 10,200 GiB free on scratch (20 tasks × 500 calls × 1 GiB per retained case, plus 200 GiB); each evaluator retains the recovery stage’s 350 GiB per-case storage floor. Slurm runs `0-19%20`. A newly completed search case with provisional loss below 13 stops further optimization only after its case, latest-completed and best-so-far receipts are written. It then runs the same fresh native selected-point and exact-repeat gates; crossing 13 is never acceptance. The target remains a goal, not a gate.

Run these commands in order when ready. `stage_torch.sh` builds the derivative, verifies the source inventory, derives and retrieves the 20-start plan, and performs zero model solves. It does not submit jobs.

```bash
bash code/cluster/estate_birth_calibration/binary_continuation_v2/stage_torch.sh
bash code/cluster/estate_birth_calibration/binary_continuation_v2/smoke_torch.sh 0
```

After the smoke task finishes, inspect its terminal receipt and run the remote receipt writer once:

```bash
ssh torch '/share/apps/anaconda3/2025.06/bin/python /scratch/td2248/projects/estate_birth_binary_continuation_20261004_v2/record_smoke_gate.py --smoke-dir /scratch/td2248/projects/estate_birth_binary_continuation_20261004_v2/results/smoke_binary_chain_0 --job-id JOB_ID'
```

The writer accepts only exactly two objective calls, a fresh passing native selected-point verification, the 14 target rows, 31 parameter rows and the exact repeat with all 17 standard plot hashes. Production submission is blocked unless that receipt passes, and a durable lock plus receipt prevents a second submission attempt. An ambiguous `sbatch` outcome leaves a `submitting` receipt for queue inspection; the helper never retries.

```bash
bash code/cluster/estate_birth_calibration/binary_continuation_v2/submit_torch.sh
python3 code/cluster/estate_birth_calibration/binary_continuation_v2/collect_torch.py --stage STAGED_REMOTE_DIR --out COLLECTION.json
```

The collector is read-only and retains incomplete task statuses. A selected result reports the complete native target-fit and parameter packets; a provisional optimizer checkpoint does not enter the accepted result table. Do not edit the existing historical `continuation/` package. No smoke or production job is submitted by building or checking these files.
