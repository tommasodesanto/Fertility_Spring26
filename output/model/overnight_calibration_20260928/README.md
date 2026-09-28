# Overnight calibration, September 27–28

## Live state

Torch job **18687184** is submitted for six gated numerical smokes. Search is
not yet approved. New stage: `/scratch/td2248/projects/fertility_night_calibration_20260928_v1`.
Nine focused controller tests pass (job18687181), including23 owned synthetic
subprocesses covering smokes, search, hourly rendering and eight final repeats.
Zero-solve source/contract verification passes. Contract/source pins and launch
receipt are in `contract_v1/` and `cluster/launch.json`.

## Authorized design

Torch-only calibration and failure diagnosis, no economic model changes.
Hard end **08:00 EDT September28** (1790596800); search ends07:00 and repeats
07:50. Queue/setup consume this fixed window. Hourly heartbeat is active.
Up to24single-thread model workers,192GB allocation; no heavy Mac work.
Maximum720search+6smokes+8finalrepeats+12reserveddiagnostics=746objectives.
Each objective has1800seconds, retaining23SSlimit and all scientific gates.
No extra diagnostic objectives are automatically launched.

Model, grids, targets, bounds, closure and normalization settings are identical
to the approved evening contract. Primary/block/relative-identity weight lanes
stay separate; rescore every successful point under common primary weights.
All three validation rows remain visible. Final selected candidates need two
fresh repeats, full14target/31parameter tables and17standard diagnostic plots.

Search starts at the repeated evening block0347 candidate. Replace unsuccessful
full-box draws with local/moderate joint, subspace and coordinate proposals;
retain full bounds. First six search cases replay two exact evening timeout
points across the three lanes under the longer wall-clock allowance.

## Failure evidence and interventions

`failure_diagnosis/` retains all360evening case classifications and stage evidence.
All77housing failures occurred in84broad draws; seven others timed out. All61
timeouts had4–8certified equilibrium solves; seven finished fertility normalization
before timeout. This does not prove equilibrium nonexistence. Lead doubled only
objective wall-clock allowance900→1800seconds; no scientific gate was relaxed.

Initial source preparation was temporarily blocked by session permissions;
author restored access. Earlier automation/SSH/Git denials remain in status history.
Preflight18687097 and diagnostic18687169 failed because the synthetic hourly
render fixture omitted its runtime reference. Test-only repair added the missing
pin; controller and household runtime unchanged. Both logs and the remote failed
fixture are retained;18687181 passes. Independent review also strengthened checks
of all saved anchor/repeat CSV hashes before contract use.

Next: inspect all six numerical smoke comparisons and17plots before explicit
lead approval releases search. Old evening run18672459 remains immutable.
