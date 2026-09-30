# Fixed-credit runtime validation v3

Reference: **2007 stationary reference — block0506, September 28 verified export**.

## Why v3 exists

The v2 smoke job **18845938** failed before Python started. Slurm copied the submitted `launch_smoke.sh` into its spool directory; line 11 used `dirname "$0"` and therefore looked for `launch_runtime_validation.sh` beside the spool copy instead of in the staged source. v2 remains intact as the failure record.

v3 copies the v2 driver and economic plan. Exact changes from v2 are: (1) both Slurm wrappers replace `exec "$(dirname "$0")/launch_runtime_validation.sh" MODE` with `exec "/scratch/td2248/projects/fixed_reference_credit_runtime_20260930_v3/source/launch_runtime_validation.sh" MODE`; (2) `launch_runtime_validation.sh` uses the v3 root and explicit v3 `/source` path; (3) `stage_source.sh` stages beneath the v3 root; (4) `plan.json` retains its v2 schema and changes only `remote_root` to the v3 root; and (5) Slurm job names identify v3. `run_runtime_validation.py` is byte-identical to v2. Time and memory budgets, driver logic, mode behavior, and economic inputs are unchanged.

## Staging and launch

**September 29 23:45 New York monitor:** smoke **18846467 COMPLETED 0:0**
in 41 seconds, zero lifecycle evaluations. Lead authenticated receipt pins,
module origins, compiled fixtures, and actual-checkpoint cash arithmetic.
The two entrant cells have cash before rent/consumption of -0.1340497885 and
-0.0677697562; buying resources are also negative. Strict-zero infeasibility
is confirmed without an entry change. Compact receipts and the lead-review
record are in `collected/smoke/` and `collected/SMOKE_LEAD_REVIEWED`.

**Authorized control successor submitted once: 18849552.** Remote dispatch
record: `results/control_submitted.txt`; log: `results/control_18849552.log`.
Do not submit another control. Read `results/control/launch.json` for its
unchangeable clocks. This replay leaves the scalar unset and performs at most
one lifecycle evaluation at the reference price; it is not a new GE solution.

**Submitted once:** Torch smoke **18846467**, September 29 evening. Immutable
sources were staged and made read-only before submission. Driver SHA256
`f9940add31491a6132cff52867d469f91b6134b1e5572b42310aea5a80caf4c1`;
plan SHA256 `28ac8cc51c05198789a782b5af18dd59d0a73b36ec89cf9c83115409da585430`.
The existing hourly heartbeat `check-frozen-reference-borrowing-ge` tracks this
version and is authorized to submit the one reviewed control successor.
Record any successor job ID here before ending that dispatch turn.
Stop after verified delivery, substantive failure or by September 30 09:00
New York. No further repair-and-retry is authorized by this packet.

Run `bash stage_source.sh` once on Torch, then submit `sbatch launch_smoke.sh` once. The staged source is immutable at:

`/scratch/td2248/projects/fixed_reference_credit_runtime_20260930_v3/source`

Results and scheduler logs belong under the sibling `results/` directory. After a passing smoke receipt has been reviewed, the authorized successor is one exact baseline control. Create `results/SMOKE_LEAD_REVIEWED` only after review, then submit `sbatch launch_control.sh` once. The control retains the v2 fixed clocks (300 seconds for the case, 900 seconds total), one thread and 24 GiB. No retries or other model runs are authorized by this packet.

Driver and economics remain byte-identical to v2. The smoke permits zero lifecycle evaluations; control performs its single exact flag-off replay only after review. The reference, overlay source pins, strict-zero affordability proof, and other gates are unchanged.
