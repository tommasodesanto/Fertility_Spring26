# Native fixed-price financing diagnostic

September 20 onward: the author has authorized a new bounded
[specification follow-up](specification_followup/README.md), with an active goal
and check-ins. Its launch manifest and reviewed findings will live there; the
completed overnight battery below remains a separate fixed experiment.

This is a full-lifecycle partial-equilibrium diagnostic from the frozen September-14 native checkpoint. It holds checkpoint/schema inspection and one receipt plus compressed native arrays for every completed arm. Every arm applies the native sequential-fertility and current-choice period evaluator to the saved beginning-of-period mass. It does not calculate a stationary endpoint, a transition, a market-clearing price, or a fiscal-clearing transfer.

Run two independent `baseline` cases first. Both require all saved policy arrays to reproduce at `atol=1e-10, rtol=0`; only then run `mortgage_only`, `unsecured_only`, and `both` as independently timed processes.

## Reviewed overnight result — September 20 morning

All final-chain jobs completed: search 18049121 (2h34m35s), original grid
18050123 (41m36s), pilot smoke/grid 18050440/18050441 (3m19s/52m), and refit
smoke/grid 18050443/18050451 (3m41s/53m23s), all exit 0. No further jobs were
launched. [Consolidated review receipt](overnight/final_review.json).

The classic persistent-plus-transitory income search selected case 60 from
96 proposals, 89 valid and seven rejected. Its loss is 353.6588729140903,
versus 502.74561411262744 for the prior local search and 179.2984242480252 for
the retained original benchmark. All 13 target definitions, values and weights
match the original comparison. This is a finite local diagnostic; identification,
convergence and adoption are not established. The housing fit remains weak:
ownership among heads 30–55 is 48.43% versus 64.83%, and capped mean rooms are
6.677 versus 5.561. The recent-parent ownership gap improves greatly, but mean
rooms worsen relative to the prior selected point. Annual discounting 0.99 and
first-child housing requirement 2.3 are at their active upper bounds.

[Full 13-row fit and 17-row parameter/bounds table](overnight/final_search/readout.md),
[source/repetition/plot receipt](overnight/final_search/receipt.json), and
[selected 17 standard plots](overnight/final_search/selected_standard_diagnostics/)
are collected. The selected normalized child-preference level is 0.1922702943;
the summary's 0.2394995040 is the initial normalization seed. Raw scorer CSVs retain
the generic beta bound 0.9995; the readout reports the actual 0.99 restriction.
Two native repetitions agree exactly on numeric fit, loss and price. Refit
mechanisms use the second verified checkpoint; its file hash differs from the
first repetition while the verified economic results agree. Lead visually
inspected fertility by age, ownership by age and market clearing; all 17 selected
plot hashes were checked, without claiming visual review of every plot.

The mechanism battery has 144 distinct cells and six extra controls. All 150
case receipts passed the original gates, and the graph manifest records 2,550
primary PNGs. Prices, preferences and initial population stay fixed within each
family; their values can differ across families. Mortgage relaxation changes
both deposit and collateral access, with financed share 0.8 to 1.

| Family | Mortgage: birth-flow change | Mortgage: explicit lifetime cohort births |
|---|---:|---:|
| Original | +0.690% | +0.234% |
| New-income pilot | +0.174% | +0.085% |
| New-income refit | +0.230% | +0.081% |

These are different outcome definitions. Birth flow uses a fixed initial
cross-section; explicit lifetime births sum all births in a simulated entry
cohort. Neither is the separate stationary 2.1 normalization. At rental cap 10,
the refit's mortgage birth-flow effect is only +0.00108%, compared with +0.230%
at cap 6. This points to a rental-access interaction, not an identified mediation
decomposition or a general-equilibrium result.

The unsecured-credit result is not robust across income specifications or outcome
clocks. With the original earnings process, the largest credit dose raises cohort
births about 3.02%. In the new-income refit, credit of 0.25 times four-year earnings
(one annual earnings amount) raises current birth flow 2.456% but lowers explicit
cohort births 0.769%; credit of one four-year earnings amount gives +2.380% flow
and -2.415% cohort births. Five four-year earnings amounts give identical new-income
outcomes to one; the reason for that plateau has not been established. Credit is
zero in retirement. These are diagnostic doses, not approved policy defaults.

A bounded saved-output check found exactly identical entry distributions across
the refit's baseline and credit arms. Summing the 17 age records exactly reproduces
lifetime births. Higher births at age 18 are outweighed by lower later births;
lifetime first births also decline. Thus the discrepancy is not an entry-population
change or a summation error. The economic or numerical cause is still outstanding.
Do not treat positive initial birth flows as proof of higher lifetime fertility.
[Accounting evidence](overnight/final_mechanisms/cohort_accounting_check.json).

[Complete contrasts](overnight/final_mechanisms/comparison.csv),
[mechanism readout](overnight/final_mechanisms/comparison.md),
[all family receipts](overnight/final_mechanisms/receipt.json), and
[full graph manifest](overnight/final_mechanisms/graph_manifest.json) are collected.
The next substantive decision is to diagnose the new-income housing fit and the
credit/cohort response before adopting the earnings change or strengthening the
slides' mortgage mechanism. The overnight battery alone supports neither a
converged replacement calibration nor a literal frictionless benchmark.

## Completion receipt — September 20, 00:02 EDT

Torch access is restored. Original production 18050123 completed in 41m36s,
exit 0. The collected [summary](overnight/monitor_20260920_0002/original/summary.json)
and [comparison table](overnight/monitor_20260920_0002/original/comparisons.csv)
contain 50 cases: two controls plus all 48 distinct combinations. Lead review
checked all case/cohort statuses, the budget and occupied-value gates, matching
baseline metrics, and exact numerical equality between comparison rows and case
receipts. Each arm records 17 primary standard diagnostic plots. The
[packet manifest](overnight/monitor_20260920_0002/graph_manifest.json) confirms
17 PNGs in each of 50 remote primary-packet directories and their total sizes;
PNGs were not recopied or newly inspected, and this is not a file-hash audit.

New-income smoke 18050440 completed in 3m19s, exit 0; its three
[case receipts](overnight/monitor_20260920_0002/stationary_new_income_smoke/summary.json)
pass the same checks and retain the pinned evaluated-population input.
New-income grid 18050441 and earnings search 18049121 are still running in the
scheduler snapshot. Refit smoke/grid 18050443/18050451 await their dependencies.
The [progress receipt](overnight/monitor_20260920_0002.json) records 35/50 pilot
cases including two controls, and 32 processed search proposals: 27 valid and
five rejected, with an active heartbeat. Rejected proposals are expected search
outcomes rather than overall job failures. No job was restarted or altered.
The final economic comparison and verified
selected search packet remain outstanding; the follow-up stays active.

## Monitoring access — September 19, 23:43 EDT

SSH authentication to `torch` failed before the scheduled remote query. Current
job states, progress counts and failure logs could not be checked. This is a
monitoring-access failure, not evidence that any submitted job failed. No jobs
were changed or restarted. All dependencies remain cluster-managed and do not
require the laptop to stay open. Restore SSH access before collecting fresh
receipts; the follow-up remains active, with no repeated notification while the
same access blocker persists. [Receipt](overnight/monitor_20260919_2343.json).

## Full battery and population correction — September 19, 23:26 EDT

| Family | Smoke | 48-cell production | Dependency |
|---|---|---|---|
| Original paper checkpoint | 18050122, passed | 18050123, running | Own smoke |
| New-income stationary pilot | 18050440 | 18050441 | Own smoke |
| New-income selected refit | 18050443 | 18050451 | Verified search 18049121, then own smoke |

All three use the dose design below: 48 cells and two controls per production
job, plus three smoke cases. The two new-income snapshots are separate immutable
folders `finance_dose_income_v2` and `finance_dose_refit_v2`, under the same
native-financing remote parent. The original `finance_dose_v1` remains untouched.
The full 144-cell mechanism battery complements up to 96 joint search proposals.

The old new-income matrix failed a bitwise input-population check. The corrected
read-only replay, job 18050183, completed in 14 seconds: raw stationary mass versus
saved `evaluation.g_pre` has L1 3.5100904979533746e-15, maximum difference
1.62824089953466e-15, 17 changed entries and no total-mass change. One saved-policy
replay reproduces the saved evaluated pre-choice mass, current mass and births
exactly. Job 18050052 was the preceding failed-import diagnostic; no model solve
ran in either diagnostic. [Evidence](overnight/population_diagnostic/report.json).

The new-income grids therefore explicitly select the saved evaluated pre-choice
mass as their common input (`--population-source saved_evaluation`). This requires
finite arrays of identical shape and absolute L1, maximum and total-mass differences
from raw stationary mass no greater than 1e-12; larger differences fail. The
per-arm `np.array_equal` population gate and original baseline/policy/budget/value
checks are unchanged. Each checkpoint records the input choice and measured
differences. The original family retains its raw stationary input. This corrects
a floating-point-scale input discrepancy; no solver, objective, preferences or
calibration result is changed. Ten focused tests, shell syntax and lead line
review passed; the candidate-family numerical smokes remain pending.

Reproduce with a fresh `DOSE_TAG`, `POPULATION_SOURCE=saved_evaluation`,
`FAMILY=stationary_new_income` or `FAMILY=refit_new_income`, and `SUBMIT=1` using
`code/cluster/submit_e5f_financing_dose.sh`; the refit requires `AFTER_JOB=18049121`.
Every job fails closed on contract or numerical failure. Each mechanism grid has
a 600-second per-case limit, five-hour controller budget and six-hour allocation.
Expected time is 50–100 minutes per family. Future follow-ups collect full search
fit/parameter tables and the unchanged diagnostic packets before interpretation.
The full battery is diagnostic; cross-family prices, preferences and native entry
wealth may differ, and it establishes neither convergence nor general equilibrium.

[All jobs and stopping rules](overnight/expanded_submission.json).

## Expanded overnight battery — September 19, 23:16 EDT

The search smoke 18049120 passed, and the 96-proposal search 18049121 is running.
The renewed author request adds a 48-cell original-family dose experiment:
financed shares 0.8, 0.9, 0.95, 1; unsecured borrowing limits 0, 0.25, 1, 5 times
four-year after-tax earnings (zero in retirement); and rental room caps 6, 8, 10.
For example, 0.25 times four-year earnings is one annual earnings amount.
The lower credit doses distinguish a response at moderate borrowing limits from
one that requires the unusually large five-period limit. Mortgage access changes
deposit and collateral constraints jointly, so these are not pure deposit effects.

Smoke 18050122 completed in 2m50s (exit 0): two exact baseline controls and the
combined extreme arm passed, each with a completed cohort and 17 standard plots.
Production 18050123 is running and covers 48 cells plus two controls.
Each case includes fixed-population total/first birth flows, housing and tenure,
explicit lifetime cohort births and first-birth age, and 17 standard plots.
Prices, preferences and initial population stay fixed; there is no GE claim or
completed-fertility normalization in this exercise. Comparisons use each family’s
own unchanged baseline. Source and numerical gates are unchanged.

The 50 production arms imply approximately 50–100 minutes at 60–120 seconds each.
Each case has a ten-minute limit, the controller a five-hour budget, and Slurm a
six-hour allocation. Progress is written every 30 seconds and after every arm;
the controller stops at the first failed gate. Seven focused tests and shell
syntax checks passed. The exact-loop cluster smoke passed.

New-income mechanism production is not resubmitted: a bounded population replay
is diagnosing the failed equality guard. The original failed matrix and cancelled
dependents remain historical receipts. The running dose snapshot is isolated at
`/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/finance_dose_v1`.
Reproduce staging/submission with `SUBMIT=1 bash code/cluster/submit_e5f_financing_dose.sh`
and a **fresh** `DOSE_TAG`; do not overwrite a staged or running snapshot.
[Expanded receipt](overnight/finance_dose_v1/submission.json) and
[pinned plan](overnight/finance_dose_v1/plan.remote.json).

## Overnight progress — September 19, 22:45 EDT

Coordinate job 18047156 completed all 16 proposals in 41m10s. Case 3 was selected
with loss 502.74561411262744 and two verified native repetitions. Full target and
parameter tables are in the [complete readout](income_search_18047156/readout.md),
with source pins and 17 standard plots in the same folder. Lead checks confirmed
all fit contributions sum to the selected loss, target/income fingerprints match,
and both native repeats agree; the fertility-by-age plot was inspected.
No convergence or adoption is implied.

Mechanism smoke 18049130 stopped after 3m43s. Both original-family controls and
the combined treatment passed. The stationary-new-income baseline failed the
saved-initial-population equality check before completing; the difference's
size and cause remain undiagnosed. Its downstream jobs 18049131 and 18049132
were cancelled without running. No automatic retry or gate relaxation occurred.
The joint earnings-search smoke 18049120 is running, independently of that
failure, and search 18049121 remains dependent on its success.
[Monitor and failure receipt](overnight/monitor_20260919_2245.json).

## Overnight experiments — September 19

The author requested more useful experiments overnight. The following jobs
are submitted; successful smoke tests are required before production starts.
The cluster owns all dependencies, so the laptop can be closed.

| Job | Purpose | Dependency | Limit |
|---|---|---|---|
| 18049120 | Two joint proposals and two selected-point repetitions | Existing coordinate poll 18047156 succeeds | 1 hour, 2 CPUs |
| 18049121 | Up to 96 joint parameter proposals and two selected-point repetitions | 18049120 succeeds | 6 hours, 8 CPUs |
| 18049130 | Exact-loop financing smoke: two controls and combined treatment in each of two families | None | 3 hours, 1 CPU |
| 18049131 | Eight financing/rental arms plus two controls in original and stationary-income families | 18049130 succeeds | 3 hours, 1 CPU |
| 18049132 | Same eight arms plus two controls at the verified overnight refit | 18049130 and 18049121 succeed | 3 hours, 1 CPU |

Search design: original structural point plus the four best distinct valid
coordinate points are five centers. Seed 20260919 produces antithetic joint
uniform directions, interleaving centers and two local scales (0.4 and 1).
Positive-coordinate log widths are 0.55 for each fertility curvature and the
bequest shift, 0.20 for housing preference, 0.35 for housing supply scale,
and 0.30 for the first-child housing requirement. Additive widths are 0.012
for annual discounting and 0.20 for the bequest level and first-birth cost.
These are search-design choices, not priors. All nine authoritative bounds
are enforced; annual discounting remains at most 0.99. The 12 scored targets,
weights, source snapshot, and separate fertility normalization remain pinned.
Known points are excluded. Search dispatch stops after 4.5 hours or 96 proposals,
or six wholly failed batches; each point has a 900-second limit, and two
repetitions reserve 30 minutes. The smoke requires a successful new proposal
and exact selected-point verification. At roughly 400 seconds per point,
96 points on eight workers plus verification is about 110 minutes. There are
at most 784 native normalization solves for production (98 evaluations times
eight); smoke adds at most 32. No optimum or adoption claim follows.

The mechanism matrix varies the financed mortgage share between 0.8 and 1,
unsecured credit between zero and five times **four-year** after-tax earnings
(zero after retirement), and the rental room cap between 6 and 10. Mortgage
access changes deposit and collateral limits jointly. Within each family,
prices, preferences and initial stationary population are fixed. The three
families are the retained paper checkpoint, the new-income stationary pilot,
and the verified refit. The middle family retains the nine original structural
parameters but has its own normalized child-preference level. Across-family
comparisons therefore need not hold preferences, prices or entry wealth fixed.
Each arm reports stationary-population birth flows, housing and tenure, and a
full native lifetime cohort. Explicit cohort births are distinct from the 2.1
normalization target. These are partial-equilibrium diagnostics, not a literal
frictionless economy, identified mediation, or GE counterfactuals.

Per-arm limit is 600 seconds, with 30-second heartbeats. Budget, mass,
probability, occupied-value, saved-population and exact baseline-control gates
remain active; all solved arms require 17 standard diagnostic plots. The matrix
has 36 solves including both baseline-family smoke loops and repeated controls;
at roughly 60–90 seconds including cohort/graphs, expect roughly 35–55 minutes
across its stages, subject to measured smoke times. Failed dependencies block
later jobs rather than triggering retries. Local verification: 17 focused tests
and both shell syntax checks passed; live numerical smoke outcomes are pending.

[Submission receipt](overnight/submission.json) and [complete staged plan](overnight/plan.remote.json)
pin job IDs, source hashes, budgets and target/parameter bounds. Raw scored
parameter tables retain the generic 0.9995 discount upper bound; readouts must
annotate the active 0.99 search restriction from the plan. No production result
has been reported yet. The existing app follow-up checks every 15 minutes and
reports meaningful completion, failure, or required action.

## Bounded earnings refit — job 18047156

Submitted and running on Torch, with exact saved-checkpoint income-payload and
zero-solve wrapper preflight passed. The first four cases have live heartbeats.
The controller tries at most 16 one-coordinate changes around the verified
651.9418197098323 pilot. It keeps all nine structural coordinates available;
all 12 scored moments and the separate 2.1 fertility normalization stay fixed.
Annual discounting remains bounded above by 0.99; the retained scorer's generic
0.9995 bound is not the search bound. Full pilot tables appear below.

This is one bounded local poll, not a converged SMM fit or adoption decision.
Search budget: 2,100 seconds; at most 900 seconds per trial; stop dispatching
with less than 500 seconds remaining. Four single-threaded workers run at most
16 trials (roughly 200–500 seconds each based on the six-solve 459-second pilot
and a closer starting normalization). The native loop allows up to eight solves
per evaluation; estimated total 36–144 solves including final verification.
The selected point is repeated twice, reserving 1,200 seconds. Hard controller
limit 3,420 seconds; Slurm limit 3,600 seconds. Contract changes abort; numerical
failures are recorded; timeout cleanup includes descendants in new sessions.
Latest completed case, best-so-far, and 30-second case heartbeats are retained.
The selected result must provide all 13 fit rows, 17 parameter rows, and the
same 17 diagnostic plots. No further search launches automatically.

Sixteen targeted tests passed, covering constructor moments, source/target
contracts, the complete batch loop, selection/repetition, budget stops,
parameter binding, and descendant cleanup. The real pilot previously exercised
the unchanged equilibrium/normalization loop; startup also preflighted the
new two-repetition contract and checked the actual checkpoint's income arrays.

Remote results: `/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/income_search_v1/results_18047156`.
Source/target hashes and bounds: [search plan](earnings_candidate/search_plan.remote.json).
Actual checkpoint audit: [receipt](income_search_18047156/income_payload_audit.json).
Launch: `SUBMIT=1 bash code/cluster/submit_e5f_income_candidate_search.sh`.
The app follow-up checks every ten minutes and reports completion or failure.

## Submission receipt (September 19, 2026)

Torch job **18034069** (`native_finance`) was submitted from
`code/cluster/submit_e5f_native_financing_diagnostic.sh` to remote root
`/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a`.
The run uses the fixed native checkpoint and exact native parameters, ψ,
pension, rebate, baseline price, and growth rate; it does not perform GE
recalibration or a transition. The five full solves are two baseline controls
followed by `mortgage_only` (φ=1, λ=0), `unsecured_only` (φ=0.8, λ=5), and
`both` (φ=1, λ=5). The first error halts the batch; there are no automatic
retries. Job **18034069** completed in **3m24s** with exit 0. Both baseline
controls reproduced saved policies, current mass, and births at
`atol=1e-10, rtol=0`; all five cases passed zero budget-excess mass and
occupied wealth-value-drop gates.

Current-calibration four-year birth flows at common initial population and
prices were baseline **0.1155903878306699**, mortgage-only
**0.11638809385408035** (**+0.690114496872%**), unsecured-only
**0.12377387987346156** (**+7.079734047419%**), and both
**0.12451625195091737** (**+7.721977828574%**). These are birth flows, not
completed fertility, and are partial-equilibrium diagnostics; they do not
measure GE offsets or revised-calibration effects. Mortgage-only relaxes both
the deposit and collateral limit. The supplementary measured report is
[report/report.md](report/report.md).

The credit cap λ=5 means five times age-specific after-tax four-year income
as constructed by `build_debt_caps`, not five annual incomes; the cap is zero
after retirement. The 82/86 taper is retained only as an old experimental
specification. The PE report uses the same baseline pre-choice mass and native
sequential mapper. The diagnostic sequence is intended to establish the
current-calibration mechanism before any author-approved structural change or
recalibration; no full modern BGM earnings specification is included.

Retained native full-fit references are
[complete target fit](../paper_baseline_sep14/replay_20260917/native_output/selected_target_fit.csv)
and [complete parameter table](../paper_baseline_sep14/replay_20260917/native_output/selected_parameters.csv).
The earnings variant and recalibration/GE closure remain outstanding; no
automatic launch follows.

## Income pilot receipts

Stationary job **18040896** completed in 8m28s with exit 0; cohort job
**18041070** completed in 2m16s with exit 0. The stationary pilot held all nine structural parameters fixed and solved
the new earnings equilibrium with normalized
\(\psi=0.2429719621740803\), yielding 2.1000283897. Its six GE solves and 17
plots passed. The objective is **651.9418197098323**, versus retained
**179.2984242480252** under the identical objective. The cohort exercise is a
different fixed-preference exercise: explicit birth flows are
**1.87241076 → 1.66643427**, first-birth probability **0.82046004 →
0.76890885**, and mean grid first-birth age **24.095 → 25.702**. Mortgage
effects are small; total births and first births are reported separately.
Do not read cohort birth flows as completed fertility or call them TFR.

The stationary target-fit table below reports all 13 rows, including the
unweighted normalization row. Displayed numbers are rounded; the linked CSVs
retain full precision. `—` denotes a null weight or loss contribution.

| Moment | Target | Model | Gap | Weight | Loss |
|---|---:|---:|---:|---:|---:|
| Initial model completed fertility | 2.1 | 2.10003 | 2.83897e-05 | — | — |
| Childless women, ages 40–44 | 0.198279 | 0.233382 | 0.0351036 | 35532.3 | 43.7851 |
| Exactly one child among mothers, ages 40–44 | 0.213655 | 0.198902 | -0.0147534 | 26952.8 | 5.86666 |
| Period mean first-birth age | 25.9763 | 27.4911 | 1.51479 | 139.828 | 320.847 |
| First births at age 30+ | 0.249278 | 0.307885 | 0.0586071 | 13866.1 | 47.6271 |
| Wealth / annual gross labor earnings | 6.14586 | 6.69313 | 0.547266 | 7.5951 | 2.27473 |
| Annual bequests / aggregate wealth | 0.0088 | 0.00700983 | -0.00179017 | 5.16529e+06 | 16.5533 |
| Old wealth/income p90 / median, ages 76–84 | 3.51594 | 4.55287 | 1.03693 | 10.6164 | 11.415 |
| Mean occupied rooms, capped at 9 | 5.5611 | 6.17848 | 0.61738 | 128.021 | 48.7962 |
| Ownership, heads 30–55 | 0.648334 | 0.490673 | -0.157661 | 2339.36 | 58.1499 |
| First-birth room response, −1 to +3 | 0.720246 | 1.29942 | 0.579172 | 137.565 | 46.145 |
| Rooms: 3+ versus 1–2 resident children (model dependent proxy) | 0.347067 | 0.350358 | 0.00329151 | 280.528 | 0.00303925 |
| Recent-parent ownership gap | 0.162896 | 0.119701 | -0.0431943 | 27055.8 | 50.4793 |

The full 17-parameter diagnostic table reports estimate, bound, near-bound
flag, and status. Beta is **0.99** at the active restricted upper bound;
the raw generic table's upper bound is **0.9995** and should not be confused
with the active pilot bound.

| Parameter | Estimate | Active bounds | Near active bound | Status |
|---|---:|---:|:---:|---|
| beta_annual | 0.99 | 0.94–0.99 | True | Retained value; fixed in this pilot |
| kappa_fert | 0.337734 | 0.02–50.0 | True | Retained value; fixed in this pilot |
| kappa_fert_continuation | 0.397856 | 0.02–50.0 | True | Retained value; fixed in this pilot |
| chi | 1.04965 | 0.1–5.0 | False | Retained value; fixed in this pilot |
| H0 | 8.1121 | 0.2–80.0 | False | Retained value; fixed in this pilot |
| theta0 | 0.081051 | 0.0–8.0 | False | Retained value; fixed in this pilot |
| theta1 | 0.0852036 | 0.02–16.0 | True | Retained value; fixed in this pilot |
| first_birth_fixed_cost | 0.265765 | 0.0–8.0 | False | Retained value; fixed in this pilot |
| h_P | 2.3 | 0.1–2.3 | True | Retained value; fixed in this pilot |
| hbar_child_rooms | 0 | fixed | False | zero restriction |
| psi_child | 0.242972 | normalized | False | normalized to 2.1 |
| payroll_tax | 0.179 | fixed | False | externally fixed |
| pension_period | 2.04636 | derived | False | budget derived |
| housing_supply_elasticity | 0.63 | fixed | False | externally fixed |
| tenure_choice_kappa | 0.005 | fixed | False | externally fixed |
| alpha_cons | 0.733 | fixed | False | externally fixed |
| sigma | 2 | fixed | False | externally fixed |

Source pins: [target-fit CSV](income_stationary_18040896/target_fit_comparison.csv),
[parameter CSV](income_stationary_18040896/parameters_comparison.csv),
[stationary receipt](income_stationary_18040896/collection_receipt.json), and
[cohort receipt](income_cohort_18041070/job_receipt.json). Bounded refit
preparation is in progress elsewhere; no full recalibration has been launched.

## Rental-access interaction and earnings candidate

Grouped native job **18037585** completed in 28 seconds and rental-access job
**18037587** completed in 3m13s, both with exit 0. At common prices,
pre-population, and fertility preferences, expanding the native rental room
cap from 6 to 10 rooms reduced the mortgage-only birth-flow increment by
**84.1755512677%** and the first-birth-flow increment by **88.4941545740%**.
This supports a rental space-access role; it does not identify mediation or
establish GE, recalibration, or a sign proof. All three rental cases passed
zero budget/mass/value violations, the baseline control reproduced exactly,
and each case wrote 17 standard PNGs. Group reconciliation was below `1e-12`.
The grouped mortgage comparison attributes 98.89% of the mortgage first-birth increment
to initial renters, who are 96.53% of eligible mass; renter Q3 first-birth
rate rose 0.412 percentage points and mean rooms rose 0.342. See
`grouped/metadata.json`, `grouped/groups.csv`, and
`rental_access/comparisons.csv`.

The persistent-plus-iid-transitory candidate uses the existing nested
13-moment, 3-parameter fit: objective `39.5085474236852` versus
`14.3178085110745` for the full fixed-effect fit. Flag: materially worse
earnings fit; no full-model recalibration has run. See
`earnings_candidate/README.md` and `earnings_candidate/candidate.json`.

## Reporting verification

Report-only job **18036306** completed in 11 seconds with exit 0. All four
comparison rows are available. First births use the exact loss of childless
mass at the fertility stage; housing uses realized current-tenure mass.
The four-panel figure is supplementary; the full standard policy-function
gallery and a separate fertility-probability normalization audit remain
outstanding. These results do not isolate housing mediation or completed
fertility.

The native solver stores tenure probabilities in float32. Their largest
occupied-state sum discrepancy is 6.69e-8, within the explicit reporting
tolerance 2e-11 plus one storage epsilon (1.1923e-7 total). Location
probabilities sum exactly to one on inspected states. This reporting check
never changes probabilities or the model's existing acceptance checks.
Earlier reporting jobs failed on precision assumptions and serialization;
none reran the household model. See `review_receipt.json`.

To regenerate the supplementary packet on Torch (Anaconda 2025.06), run the
following on a compute node, using the retained experiment and native source:

```sh
experiment=/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a
native=/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches
python "$experiment/code/model/tools/build_e5f_native_financing_report.py" --input "$experiment/run" --output "$experiment/report" --checkpoint "$native/baseline_replay_20260917/replay/case/evaluation/raw/repetition_02/initial_state.pkl.gz" --source-root "$native/final_night_20260913/corrected_initial_source_v2/code/model"
```

## Earnings cohort launch and stop

Job **18037924**, submitted with a one-hour allocation, failed after 20 seconds
at `pre-mass changed at age 0` in the baseline cohort. No candidate earnings
comparison or recalibration was produced. The feasibility projection changed
the standardized entrant distribution; correction and a cluster smoke test
are required before resubmission. See `income_cohort_18037924/status.json`.
Launcher: `code/cluster/submit_e5f_native_income_followup.sh`.

Job **18040822** replaces the independent redraw with the native conditional
entry rule for both earnings processes, following gate-only job **18040625**.
The original native entrant has zero gate-induced L1 change; the independent
redraw changes by 0.0023731346586869503. All feasibility checks remain active.
The entry asset distribution may change with income. See
`income_cohort_18040822/submission.json`; this is a submitted diagnostic,
not a completed recalibration.

Job **18040822** subsequently stopped after both old-income cohorts completed:
the first candidate solve reached the standard plot routine with an age-zero-only
distribution, causing division by zero for empty older ages. No completed
candidate comparison is claimed. The corrected plotting packet must use the
normalized full lifetime cohort cross-section.

Independent stationary pilot **18040896** passed candidate-specific frozen
wrapper preflight and is running. It keeps nine structural parameters fixed,
solves the new income specification from fresh policies/distributions, and
normalizes fertility to 2.1 under the unchanged target contract. It permits
at most eight stationary solves (1800 seconds native, 2100 seconds wrapper).
See `income_stationary_18040896/submission.json` and
`earnings_candidate/calibration_plan.remote.json`. It is not a structural search.

Latest cohort submission: **18041070**. The graph packet now uses a normalized
full lifetime cohort distribution after all ages complete, preserving the
original standard plot set. Stationary pilot **18040896** has entered its first
native equilibrium solve. Both jobs have one-hour allocations; no completed
new-income comparison or recalibration is claimed at submission.
