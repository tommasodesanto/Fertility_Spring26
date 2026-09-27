# Supervised overnight calibration — September 27

## September 27, 03:53 EDT — hourly plots checked; follow-up prepared, not launched

The 03:45 snapshot has 168 primary and 56 identity search cases, with no
failures or resource interventions. Ten local evaluators are active, observed
RSS 21.30 GiB and free disk 221 GiB. Torch login remains unavailable. The best
provisional primary loss is 87.956; all targets and fitted parameters/bounds
are in `output/model/supervised_calibration_20260927/checks/0345/`. All seventeen
standard plots for exact selected case de_0161 were inspected under
`hourly_diagnostics/0345_selected/`. Markets clear and estate funding passes;
childlessness and mean rooms improve, while early fertility stays low, mean
first-birth age moves later and wealth/earnings falls. No main-search change.

Prepared the previously announced early-fertility-weight-3000 diagnostic in
`tmp/e5f_overnight_local_20260927/portable/night_launch_v3/early_fertility_3000/`.
Its contract SHA is `859c81541054bb7ab96138d691b356c1632a21f57b87b70e038d31fd5bcc79f6`;
local authorization ID is `tommaso_authorized_20260927_local_early3000_v1`.
Only the early-fertility weight changes (100 to3000). Targets, all other weights,
parameter bounds, source, fixed inputs and normalization are unchanged. The
starting point is authenticated primary case de_0173. Two workers, at most
20 rounds, same absolute cutoffs. Six preparation and nine reporting tests
pass, as does the actual zero-solve preflight. Immutable preparer is under
`portable/followup_tools_v1/`; running tools_v4 files have not changed.

**Not launched:** wait for identity final repeats/export to complete and its
workers to exit, then run two smokes directly with the pinned driver (the
preflight folder already exists). Compare every model moment and parameter
estimate against the pinned initial case before promoting a separate contract
and starting search. Full paths, hash, environment and sizing are in
`output/model/supervised_calibration_20260927/early_weight_preparation.json`
and launch.json. Do not exceed ten concurrent local evaluators. Next plot
review due by04:45; final independent main-search repeats remain pending.

## September 27, 03:18 EDT — weight tradeoff checked; conditional next experiment

The 03:15 snapshot has 136 primary and 44 identity search cases; all pass the
retained gates, with no fatal, timeout, inadmissible or resource-stop event.
Ten local evaluators are active (23.52 GiB RSS); free disk is 229 GiB. Torch
login remains unavailable. Best provisional primary loss is 123.345; all target
and parameter tables are in `output/model/supervised_calibration_20260927/checks/0315/`.
Mean rooms improve but remain high; early fertility remains low and mean age
at first birth moves later. Final independent repeats are still pending.

A zero-solve re-ranking of saved candidates changes only the early-fertility
weight from 100 to 1000/3000/10000. Weight 1000 selects the same candidate;
weight 10000 selects early fertility 0.570 versus target 0.810. These are
re-ranked existing points, not recalibrations or evidence of an unreachable
target. CSV/JSON method and results are retained beside the fit tables.

Under the author's authorization for weight experiments, once the identity
lane has completed its final repeats/export, use its two freed workers for a
separate early-fertility-weight-3000 sensitivity, starting from the then-validated
primary best. At most twenty rounds of two workers, two exact smoke evaluations,
a fresh immutable contract/provenance fingerprint, all other primary weights
retained, and the same 07:13 search / 08:43 total cutoffs. No launch yet and no
main-search change. This moment informs child-benefit curvature, first-birth
cost and both fertility taste scales jointly with retained fertility moments;
no identifying moment is dropped. Weight 3000 is diagnostic, not an estimated
optimal weight. Last full plot review 02:45; next due by 03:45.

## September 27, 02:48 EDT — hourly economic inspection; search continues

The 02:45 table snapshot has 108 primary and 34 identity search evaluations
(excluding four acceptance cases), with no failures or resource interventions.
Ten local evaluators remain active; observed RSS 15.20 GiB and free disk 235 GiB.
Torch login remains unavailable and cluster search has not launched.

The selected provisional primary loss is 147.412. Complete targets, gaps,
weights, contributions and all fitted parameters/bounds are in
`output/model/supervised_calibration_20260927/checks/0245/`. All seventeen
standard plots were inspected for this exact selected case, de_0108, under
`hourly_diagnostics/0245_selected/`. Market residual is 3.867e-6 and estate
funding passes. Ownership rises with age; owner demand remains concentrated
in the ten-room product. Average rooms fall toward target, but remain high.
Early fertility is moving farther below target and the first-birth housing
response remains too large: lower loss does not improve every target. Keep
this weighting tradeoff explicit in the final memo. No model, weight or search
intervention; final independent repeats remain pending. Next full plot check
due by 03:45.

## September 27, 02:17 EDT — 104 search cases completed without failures

The 02:15 snapshot has 80 primary and 24 identity search evaluations, plus four
acceptance cases in the comparison collector. Ten local evaluators are active;
RSS is 18.41 GiB and disk availability 241 GiB. No fatal, timeout, inadmissible or
resource-stop event is recorded. Torch authentication remains unavailable.

Best provisional primary loss is 202.087. Complete target and parameter tables
are in `output/model/supervised_calibration_20260927/checks/0215/`. Mean first-birth
age now matches closely, but mean rooms remain high and childlessness/early
fertility low. The identity objective's own winner has curvature 0.008 and a
common-primary loss of 635.407; raw identity-loss improvement does not indicate
superior primary fit. No model or weight intervention. Last full standard-plot
inspection was 01:45; next is due by 02:45. If the fixed thirty-round budgets
finish before the 07:13 search cutoff, evaluate a new bounded continuation from
validated best points, with a fresh immutable contract and smoke; do not alter
running contracts. Independent final repeats remain pending.

## September 27, 01:48 EDT — local fit improves; no search failures

At the frozen 01:45 table snapshot, 43 primary and 12 identity search cases have
completed, excluding the four corresponding acceptance cases. No fatal error,
timeout, inadmissible result or memory intervention is recorded. Capacity remains
eight primary plus two identity workers (nine active at the health snapshot as a
batch finished); aggregate evaluator RSS was 17.55 GiB. Torch authentication is
still unavailable and no cluster search has launched.

The best provisional primary loss is 289.289, versus the starting 592.815.
Full fourteen-row fits and all ten fitted parameters/bounds are in
`output/model/supervised_calibration_20260927/checks/0145/`. Ownership and the
recent-parent ownership gap are close to target; average rooms, childlessness
and early fertility remain material misses. This selected point has not yet
received the end-of-search independent repeats. All seventeen standard plots
were inspected for nearby candidate de_0032 (loss 290.245), distinctly from the
snapshot winner de_0044 (289.289); markets clear, ownership rises with age and
owner demand remains concentrated in the ten-room product. No economic change
or search intervention is warranted by this check. Identity results remain
separate and are also scored using common primary weights.

## Current state,01:07 EDT

Ten local searchworkers launched:8primary+2identity. Six complete acceptance
evaluations passed (primary/cold/identity pairs); each repeats all14targets and
31parameter rows exactly andwrites17standardfigures. Primary startingloss592.815
and all tables are in `acceptance_start/`; this is a starting point, notfinalfit.
Warm/cold comparisons are in `local_warm_cold/`:13vs19priceevaluations,231vs312
seconds, maximum moment difference<1e-5. Allmarket/fiscal/funding/policygates pass.
The starting17diagnostic figures were visually inspected via two contact sheets.
Ownership rises with age/income; owner demand is concentrated at the largest
retained product; meanhousing andparentownership gaps are keyeconomic misses.
Keep inheritedwide-wealth-axis/overlappinglegend limitations distinct fromcode
or economic failures; stable17graphs remain unchanged.

Torch authexpired afterjobs18624270/18624271. ClustersearchNOTsubmitted;
old-floor replay andcross-host arraycomparisons awaitcollection. Localbounded
EXPLORATION proceeds underexplicit author delegation after new-specification
acceptance; no old-baseline certification is asserted. User asked torefreshSSH.

Localimmutable tools_v4/night_launch_v2 paths are inlaunch.json andportable
local_registry.json. Primary/identitysearchwrapperPIDs7730/7731. Memoryguard7737
has32GiBaggregate-owned-evaluator threshold andrecords resourceinterventions;
it never treats astop as an admissible modelresult. Firstcancelledsmokes and
generator-pin-drift failure are preserved; no live contract was silentlyrepinned.

Author instruction at 00:42 EDT authorizes an eight-hour supervised exploration,
24 Torch workers and10 local workers, numerical/controller repairs without
changing the model, and separately labeled weight experiments. This supersedes
the previous launch hold and Torch-only restriction for this run. Deliver a
two-page PDF memo by09:30 EDT September27, with full target/parameter tables,
errors, economic misses, and next steps; detailed diagnostics stay separate.

Global numerical cutoff: **08:43 EDT /12:43 UTC**, including final repeats/export.
Search reserves the last90minutes for repeats/export. Current design: cluster
20 primary +4 identity, local8 primary +2 identity, single-thread workers. Local
resource pressure is monitored; any concurrency reduction must be disclosed.
Identity means unit weights on raw moment gaps, not an optimal or unit-invariant
criterion. Compare candidates under common primary weights as well.

Ten fitted parameters comprise nine outer coordinates plus the child-benefit
level normalized to completed fertility2.1. Fourteen target rows comprise13
scored moments and that normalization. No target is dropped. Primary retains
early-fertility weight100; curvature search bounds[0,.8] are provisional lead
choices for this authorized exploration, not final scientific adoption.

The reviewed reference is `output/model/calibration_code_integration_20260927`.
Claude's second review cleared numerical smoke, confirmed the mathematics and
previous fixes, and identified the fatal-stop export risk. Controller changes
preserve strict failure classification and will label any recovered export as
an interrupted run. Six full local acceptance evaluations subsequently passed; the old-floor replay and cross-host verification remain pending collection.

Run registry: `launch.json` (written as submissions occur). Monitoring history:
`interventions.jsonl`. The heartbeat `supervise-fertility-calibration-tonight`
checks every30minutes. At least hourly inspect actual saved diagnostic plots
and complete target/parameter tables. Missing progress for30minutes triggers
diagnosis. Immutable corrected snapshots preserve every old failure receipt.

No model changes relative to the reviewed no-floor compensated first-child
specification: retained grids, income/entry distributions, credit, conception,
preferences except estimated coordinates, and all numerical gates. Existing
measurement mismatches remain visible. Runtime portability and controller
changes are numerical/workflow changes. Numerical tests, runtime sizing and
launch receipts are added as available; preparation is not a completed search.
