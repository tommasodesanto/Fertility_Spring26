# Revised-interest timing: original-target continuation

This package extends the verified post-interest transaction-timing calibration
with ten additional Nelder–Mead chains. Tommaso adopted that timing as the
working standard on October 3. This continuation keeps the original 14 target
rows and weights, including the PSID wealth-to-earnings target
`6.92658379107299`, and the original ten free-parameter bounds, including
annual `beta_annual` in `[0.94, 0.99]`. The narrower 4.458 wealth target and
the lower 0.93 beta bound belong to a separate experimental cohort.

The ten starts are the ten lowest-loss distinct, numerically verified revised-
timing endpoints in the completed 48-chain comparison. Chain 0 starts from
previous alternative chain 13 (native loss `13.771131463467462`). Every
source `completed.json`, its SHA-256, the full start vector, original soft
source anchor, target/weight fingerprints, and bounds are pinned in
[`start_plan.json`](../../../../output/model/fixed_reference_economics_20260928/soft_timing_continuation_20261003_v1/start_plan.json).
These earlier endpoints and the new search do not have optimizer convergence
certificates. A lower loss is a numerical result, not proof of identification.

The source archive is built from the passed October 2 v3 archive, plus the
continuation driver, start plan, and ten verified endpoint receipts. No dirty
working tree upload, economic fallback, target edit, or source-model edit is
used. Each chain has one CPU, 24 GiB, six hours, at most 500 objective calls,
and a final 1,800-second native validation reserve. The driver writes
`heartbeat.json`, `latest_completed.json`, and `best_so_far.json` at every
completed case, then checks the selected point in a fresh interpreter with a
14-row fit, 31-parameter record, exact repeat, and 17 standard plots.

The isolated Torch root is
`/scratch/td2248/projects/soft_timing_continuation_20261003_v1`.
Stage and preflight from the repository root with
`bash code/cluster/soft_timing_calibration/continuation/stage_torch.sh`.
Stage verification is zero-solve. Run one chain-0 mock loop and one chain-0
two-evaluation native smoke as separate Slurm jobs:

```bash
ssh torch 'cd /scratch/td2248/projects/soft_timing_continuation_20261003_v1 && sbatch --parsable --array=0 --time=01:30:00 --export=ALL,CONTINUE_RUN_MODE=mock launch_torch.sh'
ssh torch 'cd /scratch/td2248/projects/soft_timing_continuation_20261003_v1 && sbatch --parsable --array=0 --time=01:30:00 --export=ALL,CONTINUE_RUN_MODE=smoke launch_torch.sh'
```

`verify_smoke_gate.py` requires the full native smoke on this exact stage.
Production is a separate lead-reviewed action: run `submit_torch.sh` **once**
on Torch after the smoke gate passes. It submits `0-9%10` with the budgets
above and refuses a duplicate submission receipt. The collector can fetch a
completed cohort with `python3 code/cluster/soft_timing_calibration/continuation/collect_torch.py
--mode production --root output/model/fixed_reference_economics_20260928/soft_timing_continuation_20261003_v1/collection
--fetch`; it requires all ten terminal native packets and rejects mixed target
or weight fingerprints. It neither retries nor adopts an estimate.
