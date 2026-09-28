# Overnight calibration, September 27–28

## September28,00:13EDT — early-fertility tradeoff diagnostic queued

Author explicitly requested focused runs to measure early-fertility improvement
and sacrifices elsewhere. Torch18690757 is PENDING, explicit StartTime07:00EDT,
EndTime07:45,10CPU/64GB. Eleven controller tests pass18690655, including actual
ten-process batches and injected anchor mismatch rejection (zero model solves).
Two exact14/31 anchor replays plus eight stronger fertility-scale proposals,
including two labeled zero-cost/loading boundary diagnostics. No model change,
no adoption, original bounds/normalization/scientific gates. Other-target loss
excludes early fertility's own contribution. Full tables/17plots retained.
Main search unchanged. Max main8+frontier10+price1=19 concurrent workers after
07:00. Ten frontier plus two price diagnostics fill12reserved slots; total746
unchanged. Per-case1800s, total2400s, hard07:45; no retries/extensions.
Plan SHA7c6db80402eb34a6368b0c3d8b19c344c4c9393fa020fa1f5fc3574262781248.
Evidence: early_frontier/README.md, approved_plan.json, tests_18690655.log,
launch.json. Diagnostic only, not calibrated SMM or global reachability test.
Hourly automation updated; final08:00 report includes this tradeoff.


## Latest check:23:55EDT

48/48 completed searches pass;24workers active; peak68.2GiB/192GB. Both old
timeout points now pass in all3lanes with longer time allowance. Full current
fits/parameters, primary-weight contributions and review:`cluster/check_2351/`.
No intervention. Final repeats remain pending. Latest available17smoke plots
reviewed again; newbest hourly exports not yet due.

## Live state

Torch job **18687184** is running search with **24 active single-thread workers**
oncs740. Six numerical smokes pass: complete14target/31parameter tables match
both their lane repeat and the saved evening block0347 physical results exactly.
Lead inspected all17standard plots; all six packets are byte-identical. The real
saved-render backend also passes (18687549), regenerating the same17PNGs and
unchanged14/31 tables with zero equilibrium solves. Lead approval SHA:
`38a2b7fbcada0759a671d1048e0be7b50dd8c02545c250d642ef7c2eb2f9e316`.

New stage: `/scratch/td2248/projects/fertility_night_calibration_20260928_v1`.
Nine controller tests+zero-solve verification pass18687181. Receipts, full smoke
tables and plots are in `cluster/smoke_review/`; all launch pins in`cluster/launch.json`.

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

Old evening run18672459 remains immutable. Starting common-primary loss is
26.368192413904882. Starting block objective2.5279564454651906 and identity
0.15368278194839702 are different objective scales, never directly compared.

## Scheduled failure diagnosis

Job18687833 is queued explicitly for September28 **07:00EDT**, ending no later
than **07:45EDT**. It waits until main search dispatch has ended and confirms
at most8main model workers before launching one diagnostic worker. Two sequential
GE calls count toward2of12reserved diagnostic slots; no new search allowance.
Same failed evening0041block point, same fixed child benefit and allmodelinputs;
only original versus half starting price differs. It records returned residuals
before unchanged housing/PAYGO gates. Results are diagnostics, not calibrated
solutions. No fertility normalization or policy certification is claimed.

Ten pure process/pin tests pass18687536. Lead reviewed the call against the first
native normalization GE. Approved immutableplan and source pins are in
`failure_diagnosis/approved_plan.json` and`approved_manifest.sha256`; remote files
under`<stage>/price_diagnostic_v1`, outputs`run_v1`. Caps1200s/arm,2450sparent,
07:45globalend. Pending18687806 was cancelled before execution because Slurm
interpreted a unitless relative begin interval as seconds; replacement18687833
uses explicit07:00EDT. No diagnostic solve was launched by that mistake.
