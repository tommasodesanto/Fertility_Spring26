# Supervised overnight calibration — September 27

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
