# Four demographic cases: overnight rerun plan

**Author decision, September 12–13, 2026. Planning record only; no launch from this request.**
Run both demographic mechanisms with and without migration. Zero migration is
the author's intended forecast specification; migration-on cases are controlled
comparisons, not automatic fallback baselines. Within the zero-migration cases,
prefer surviving-maturation entry if verified; retain the current demographic
mechanism as the fallback, with its missing orphan-care treatment disclosed.
Completion or a better fit alone does not select a mechanism or authorize
restoring migration.

## The two demographic mechanisms

| Branch | Demographic treatment | Presentation role |
|---|---|---|
| A: retained demographic treatment | Keep today's historical head-age conditioning, 2023 person anchor and subsequent person-cohort/headship law. The household operator removes dependents when the parental household exits, while the separate population calculation preserves its own surviving child cohorts. No explicit consumption or housing expenditure is assigned to these orphans. | Fallback, with this omission disclosed; do not describe it as an implemented or financed care system. |
| B: entry from surviving maturation | Dependents die when their parental household exits. Apply the existing stochastic maturation process only to surviving dependents. Domestic household entry comes only from that surviving maturation flow, with an explicit conversion/formation rule and separately counted migration. | Preferred candidate, conditional on numerical and empirical checks. Joint family mortality is an assumption, not a demonstrated empirical improvement. |

Keep source snapshots, output folders, target receipts and checkpoint chains
separate. Never splice A's initial state, demographic anchor or terminal
equilibrium into B's history. In B, the separate birth queue and person/head
rescaling cannot silently recreate entrants or overwrite the new closure.

## Four cases: separate demographic mechanism from migration

**Author addition, September 13:** retain the migration comparison while making
zero migration the intended baseline. The saved September 12–13 A continuation
has positive, externally supplied net migration in every post-2023 window;
it is not evidence for a closed-population path.

| Case | Demographic mechanism | Migration in the post-2023 forecast | Role |
|---|---|---|---|
| A0 | A: separate person cohorts and headship | Zero in every age-sex cell and year | Closed-population fallback |
| A+ | A: separate person cohorts and headship | Retain the supplied dated net-migration path | Controlled migration comparison |
| B0 | B: entry from surviving maturation | Zero outside migration/entry | Preferred closed-population candidate |
| B+ | B: entry from surviving maturation | Same external person-migration path, with an explicit conversion into B's states | Controlled migration comparison |

The migration split applies from the 2023 forecast origin onward, including
that future segment of forecasts formed before 2023. Keep the historical rule
fixed within each migration pair. A retains its observed historical household
conditioning through 2023 and its person anchor; do not describe that imposed
history as an endogenous zero-migration prediction. B must preserve its own
surviving-maturation rule rather than import A's reweighting. Each case refits
the successive-surprise history because changed future demography can change
earlier choices.

In A0, keep births, age-specific survival and the declared fixed headship
profiles, but set the person-law net-migration arrays exactly to zero. In B0,
set migration and every outside-origin entry flow exactly to zero; domestic
entry must come only from surviving maturation. Neither case may use a hidden
outside-entry residual to offset low fertility or to stabilize population.
A0 still imposes headship-based formation and dissolution, so zero migration
does not make its dependent and person records fully linked.

In the migration-on comparisons, the input is a dated net flow (arrivals minus
departures), not a constant immigration fraction or a residual refitted to the
model's births. For A+, preserve its existing empirical-source construction.
B+ needs an explicit age and unit conversion from this same person input into
household/dependent states; assigning all migrants to youngest entrants or
using an arbitrary constant is not an approved shortcut. This mapping is an
outstanding preflight item, not an implemented migration mechanism.

Where initial equilibrium contracts are identical, A0/A+ share one verified A
initial calibration and B0/B+ share one verified B initial calibration. Compare
migration switches at fixed initial states, structural parameters and preference
paths as well as after refitting shocks. Give these fixed-preference diagnostics
a separate bounded budget; do not interpret a comparison of differently fitted
histories as the isolated effect of migration.

### Closed-population terminal condition

Zero migration is not a one-line production change. The retained person-law
terminal solver uses fixed positive migration to support its stationary scale.
With no migration, subreplacement reproduction need not admit a positive
stationary population. Do not reuse an open-population terminal equilibrium,
restore migration, or reset household mass to obtain a positive endpoint.

Closed cases may follow shrinking finite-horizon paths. They still need an
explicit continuation-value boundary consistent with the zero-migration
assumption and horizon checks on the dates being reported. A positive closed
steady state is not a prerequisite for studying such a transition, but it
must be established before any stationary comparison is claimed. The boundary
construction and its fiscal consistency are outstanding implementation items;
100 dates alone do not certify them.

## Shared economic specification

- Use the current sequential-choice model, revised utility, stochastic child
  maturation and nine-parameter structural search. Do not reopen nesting or
  utility experiments tonight.
- Estimate annual beta with its existing cap of 0.99; do not fix it at the cap.
  Preserve the remaining approved bounds and external restrictions.
- **Rebate all property-tax revenue equally in every case.** Re-solve the
  initial equilibrium, every forecast and its terminal continuation with this
  same fiscal rule. PAYGO pensions balance separately at each date.
- Refit the approximate pre-2007 stationary economy using the existing twelve
  scored empirical moments, their definitions and weights. Keep the first-birth
  room response at 0.7202462623815278. No target dropping or reweighting.
  The author additionally requested a separate experimental initial calibration
  adding fertility stocks by age; see section2a. Its expanded objective must
  have its own target-and-weight fingerprint.
- Fit the four fertility windows from 2007 through 2023 using successive
  unexpected preference changes. At each change, households expect that level
  to persist; carry the realized household state forward. Hold preferences
  constant after 2023. These are fitted deterministic surprises, not an
  estimated stochastic preference process.
- The complete 2023 Data/Model table remains an untargeted validation exercise.
  Preserve its empirical-vintage labels and all thirteen moment families.

## Launch sequence and parallel work

### 1. Exact-loop smokes for all four cases, before large searches

For A, verify the equal-rebate initial, terminal-boundary and dated fiscal
calculations. For B, first exercise the demographic operator at saved policies,
then connect it to the same fiscal and equilibrium solver. No additional
household states are required for the aggregate domestic-renewal experiment;
that does not certify B+'s separate migrant allocation.

Exercise both migration settings with the exact 6-, 24- and 100-period loop
structures. The 24-period case must produce its own fitted history and policy
packet, not just a numerical starting guess for the 100-period case.
For A0/B0, every migration and outside-origin entry field from 2023 onward must
be exactly zero; an aggregate zero concealing offsetting age-cell flows is
insufficient. Audit earlier historical conditioning separately.
For A+/B+, pin the supplied flow and its units and verify the mapped entries.
Report births, child survival/dependency, maturation, household formation,
migration, population growth and fiscal flows separately. Complete one bounded
migration-switch comparison before multiplying full historical searches.

For B, with post-birth dependent stock C+, dependent deaths D, surviving
maturations M, household exits DH and household entry E, verify
`C_next = C+ - D - M` and `H_next = H - DH + E`, adding explicitly recorded
migration where applicable. Check nonnegative mass, stochastic probabilities,
terminal exits, birth units, no duplicated entry and several successive steps.

**Before B's long run:** pin how literal dependent units, the representative
3+ birth bin, mature persons and household heads map to one another. Keep
conversion, retention and migration visible as estimated, empirically
normalized, externally fixed or outstanding. Do not turn the diagnostic
outside-entry share or a mechanically calculated conversion into an empirical
production input.

Keep the author-requested initial completed-fertility normalization at 2.1
in both mechanisms. This does not itself certify replacement under either
mechanism's units. For any initial economy claimed to be closed and stationary,
verify actual reproduction: domestic household entry must equal household exits.
B's mortality, maturation and formation restriction must reconcile that check
with its fertility normalization; do not force both by undocumented rescaling
or by silently changing the target. Keep the birth-per-household versus
birth-per-woman mapping visible as a separate measurement limitation. Its twelve
scored empirical rows remain unchanged. An unresolved B gate must be reported
without preventing independent work on an admissible A case.

Every exact loop must write checkpoints, full observations, failed-case
receipts and the standard diagnostic packet before its search is admitted.

### 2. Full initial calibration for A and B, in parallel

Use all nine estimated coordinates and the complete objective. Proposed search
ceiling per branch: eighteen plus/minus coordinate probes, then up to two
eighteen-candidate waves of joint proposals informed by the complete objective
and local response matrix. These joint waves may move all nine parameters.
Keep the best valid point and perform two independent exact repetitions.
No branch wins from a favorable subset of moments.

Retain two initial calibrations if the A0/A+ and B0/B+ initial contracts are
identical: the migration experiment begins in their forecasts. Verify this
rather than assuming it from equal starting parameters. If migration must also
change an initial economy, record the distinct initial contract and revise the
budget before launch; do not silently add two structural searches.

Use up to eighteen concurrent calibration workers across all cases in total,
subject to actual Torch memory/account limits; each has one numerical thread.
Prioritize zero-migration feasibility and initial calibrations. Start a
historical pipeline from the first verified rebated initial candidate in each
branch while the structural search continues. Freeze its source and parameters.
A better initial calibration can start a separate history; it cannot silently
replace the starting distribution of an existing fitted path.

### 2a. Extra experiment: improve the inherited fertility profile

**Author addition, September13:** supplement the initial calibration with
pre-2007 fertility stocks by age. The purpose is to start younger and older
cohorts with more realistic accumulated births before fitting subsequent
preference shocks. This is an extra calibration specification, not an automatic
replacement of the original objective or a direct reweighting of saved cohorts.

Use pooled June2004/2006 CPS observations, consistent with the existing initial
fertility source. Proposed additional rows are mean children ever born at
ages25–29,30–34,35–39,40–44 and childlessness at25–29 and35–39. Inspect the full
0/1/2/3+ distribution in all these groups, with ages20–24 as a supplemental
boundary diagnostic. Existing childlessness and exactly-one-among-mothers at
40–44 remain in the original objective; do not duplicate them as new rows or
count all shares summing to one as independent information.

These observations help discipline first-birth cost and dispersion, later-birth
dispersion and the other utility/housing parameters through their implications
for birth timing and family size. They do not separately identify every
parameter or guarantee that a stationary age profile can reproduce actual
cohorts with different pre-2007 histories. Completed fertility at40–44 remains
distinct from final lifetime fertility and the replacement normalization.

Before this experimental search: reproduce the existing40–44 empirical rows;
record the source, sample, survey weights, years, age bands, uncertainty and
new objective weights; verify the stationary model age projection; and
reconcile literal child counts with the maintained3+ representative. Preserve
raw/capped empirical means and model-coded means separately. A synthetic
weight must be identified as synthetic, not presented as an empirical standard
error. Correlated/overlapping stock moments require explicit treatment.

Keep all twelve original scored rows and their weights. Score both the original
and augmented objectives on the same candidate, and report the complete tables
and parameter bounds. Re-estimate structural parameters to improve the initial
profile; do not merely alter the distribution while calling it a stationary
equilibrium. The2023completed-fertility path remains an outcome to assess
whether this improved initialization actually helps.

Keep this addition bounded: at most six full candidate evaluations per
demographic branch, including its starting point, and two selected-point
repetitions if a candidate is usable. This is at most sixteen additional
single-repetition evaluations across A and B, within the existing initial
search time budget. The migration split does not by itself duplicate this
experiment: reuse an augmented initial calibration across its two migration
cases only when their initial contracts match. The twelve primary
historical/horizon tracks remain the priority. A better augmented initial
calibration can seed a separately pinned history within the existing chain
budget; do not automatically double all
long-horizon jobs or splice parameters into a history already fitted.

Candidate empirical extraction and exact reproduction code belong under
`output/model/e5f_matched_pf_20260909a/design_research/fertility_contract/age_profile/`.
Point extraction alone does not activate a scored target or certify a model fit.

### 3. Historical shock fitting and horizon checks

**Author refinement, September 13: retain 6-, 24- and 100-period forecasts for
each of the four demographic/migration cases.** The intermediate horizon now
has its own fitted history and policy deliverable. Four economic cases therefore
have twelve historical/policy tracks; horizon length is a numerical check, not
a fifth demographic mechanism. The extra age-profile calibration in section2a
is a separate bounded initial-calibration experiment, not an automatic doubling
of these twelve tracks.

| Historical/policy track | Demographic/migration case | Explicit forecast |
|---|---|---|
| A0-short | A, zero migration | 6 four-year dates, 24 years |
| A0-intermediate | A, zero migration | 24 four-year dates, 96 years |
| A0-long | A, zero migration | 100 four-year dates, 400 years |
| A+-short | A, supplied migration | 6 four-year dates, 24 years |
| A+-intermediate | A, supplied migration | 24 four-year dates, 96 years |
| A+-long | A, supplied migration | 100 four-year dates, 400 years |
| B0-short | B, zero migration | 6 four-year dates, 24 years |
| B0-intermediate | B, zero migration | 24 four-year dates, 96 years |
| B0-long | B, zero migration | 100 four-year dates, 400 years |
| B+-short | B, supplied migration | 6 four-year dates, 24 years |
| B+-intermediate | B, supplied migration | 24 four-year dates, 96 years |
| B+-long | B, supplied migration | 100 four-year dates, 400 years |

These are numerical forecast horizons with terminal continuation values, not
different assumptions about how long the preference shock lasts. Preferences
are expected to remain at the latest revealed level in every track. There are
two initial calibrations when the migration-paired initial contracts match.
Each demographic mechanism's six forecast tracks then begin from the same
verified initial candidate. All twelve tracks fit their shock sequences
separately and write separate policy results.

Start each of the twelve tracks after its own required smoke tests, within the
shared worker and time limits. Long tracks need not wait for shorter fits;
short and intermediate results do not wait for long convergence. Reuse valid
shorter-horizon price and fiscal paths as numerical starting guesses when
available, but re-solve and refit under each track's own horizon contract. Allow
at most one additional bounded numerical-start attempt per track when useful.
Within each history, windows must run sequentially because
the next window inherits the preceding realized state. Use bracketed searches
over preference levels, saved price/pension/transfer guesses and exact replay.
Reject a failed trial and continue with the next admissible proposal; do not
terminate unrelated chains or admit a failed inherited state.

Short tracks are deliberately provisional finite-horizon experiments. Compare
both refitted paths and fixed-shock replays from the same inherited state:
refitting fertility alone could conceal a change in fitted preferences caused
by the horizon. Report differences in fitted shocks, prices, housing and policy
effects. Twenty-four dates (96 years) provide an independent intermediate fit
and policy result and can also seed the 100-date solve. A twelve-date numerical
bridge, if needed, is separately budgeted diagnostic work, not a fourth
deliverable horizon. A stalled long solve does not block the short or
intermediate track's fitted path and policy readout. An
endpoint-distance failure or material change with the horizon stays visible;
none of the three horizon labels certifies convergence. A numerical
horizon extension must hold future demographic primitives and fiscal rules
fixed under the same stated extrapolation, rather than change the economy
along with the horizon.
Retain the existing fertility tolerance and numerical gates.

### 4. Policy computations in the four cases

The priority comparison is the baseline 1% annual property tax with equal
rebate versus a 2% annual tax with equal rebate. Compute stationary comparisons
only where that case admits a verified stationary equilibrium. Do not require
a positive closed endpoint or substitute an open endpoint for the zero-migration
finite transition. State the terminal-boundary approximation and horizon evidence
for each transition. Stationary comparisons, when available, are not effects
from the inherited 2023 economy.

As soon as a demographic/migration/horizon track has an admissible fitted 2023
state and baseline continuation, run its rebated-tax policy transition from
that same state.
Compare births, population under that branch's definition, housing services,
ownership, prices/rents, consumption and both fiscal budgets. Baseline and
policy require matching horizons, demographic closure and migration inputs.
Prioritize A0/B0 policy results; A+/B+ diagnose sensitivity to migration.

Supply +20% and dependent-child LTV95% are the next independent policy jobs if
the main historical/tax work passes and time remains. Policy failure must not
erase completed calibration, history or stationary evidence. The present
six-date, horizon-unverified path cannot be silently promoted to a production
policy benchmark.

### 5. Automatic collection and morning decision

Write one review packet per demographic/migration/horizon track (twelve primary
packets) and a comparison of the four cases at each horizon containing:

- every initial target, model moment, gap, weight and loss contribution;
- every estimated parameter, bound and bound proximity, plus external inputs;
- all four fertility windows, fitted preferences, model/data path and saved
  continuation, with clear forecast-horizon status;
- all thirteen 2023 Data/Model rows with source vintages;
- the stable policy-function, lifecycle, intergenerational-allocation, price,
  quantity, market, population and fiscal diagnostic graphs;
- policy differences from a matched baseline, separating stationary evidence
  from transition evidence, and a concise list of failed or outstanding gates.

Prefer B0 only if its accounting, initial reproduction where claimed, household
choices, fiscal and market checks pass and its fit and horizon evidence support
the stated claims. Otherwise use A0 only to the extent its own checks permit,
explicitly stating the orphan-care and dependent/person-record omissions.
A+/B+ remain migration sensitivities even if they fit better or finish first.
If neither zero-migration case passes a transition/horizon gate, report that
limitation; do not substitute a migration-on case as the intended baseline.
Do not imply uniqueness or guaranteed existence.

## Compute budget and unattended operation

The proposed envelope is twelve hours from launch, with the final hour reserved
for reproduction, collection and figures. Initial search should use at most
three hours; it must not delay the first viable historical pipeline. Expanding
to four demographic/migration cases and three horizons does not extend this
shared wall-time budget or the eighteen-worker concurrency cap. Historical chains have explicit trial
and per-forecast limits and a shared deadline;
reserve the final three hours for admitted policy runs and verification. These
are caps, not forecasts of successful completion.

Initial-search ceiling, conditional on the verified reuse of two initial
calibrations: 54 objective evaluations plus two smoke repetitions and two
selected-point repetitions per A/B mechanism, or 116 single-repetition
evaluations across both. The four-case forecast-loop smokes are separate work
and must also be costed before launch. With at most eight stationary solves per
initial evaluation, this means at most 928 stationary solves. At the previously observed 2–5 minutes per
stationary solve, that is roughly 31–77 CPU-hours, or 1.7–4.3 hours with eighteen
fully utilized workers, before startup, uneven work, failures and queueing.
Replace this rough estimate with the exact-loop measured cost before submission;
reduce round counts if the three-hour stage budget requires it.

The optional age-profile experiment adds at most sixteen single-repetition
evaluations, or128stationary solves at the same eight-solve cap. Including that
pilot, the combined initial-work ceiling becomes132evaluations and1056stationary
solves (roughly35–88CPU-hours at the same rough timing). Keep the three-hour
initial-search deadline; share workers or reduce rounds rather than extend the
critical path silently.

Historical fitting has twelve primary chains (two demographic mechanisms times
two migration settings times three horizons), four windows and six preference
trials per window: at most 288 primary forecast attempts. Preserve at most one
bounded alternative numerical-start attempt per track, adding twelve attempts
for a ceiling of 300. These replace the earlier eight-track ceilings of 192/200;
they do not permit replaying whole searches. This is a 50% increase in the
attempt ceiling, not a claim of 50% longer runtime: horizon costs differ.
Do not spend the reserve on identical failed restarts. Native smokes,
fixed-shock migration/horizon comparisons and the optional augmented-initial
history each need their own recorded budget within the common deadline.
The common deadline will usually bind much earlier than these ceilings.
A six-date converged forecast previously took 6–45 minutes at fixed preferences; longer forecasts
and repeated shock roots can dominate total time. Parallelism helps independent
chains and policy branches, not the dependence between historical windows.
Before launch, record horizon-specific solve counts, memory and wall-time
estimates from the native smoke; no uncosted 100-date search. A 100-date mapping
and a 24-date mapping need separate timing receipts. A 24-date mapping has
roughly four times the household-date work of six dates: naively 10–13 minutes
per mapping using the observed short timing, before changes in iteration count.
This is an unmeasured planning estimate, not the cost of a complete fitted path.
A 100-date mapping
has about 16.7 times the household-date work of a six-date mapping. Naively
scaling the observed 2.5–3.3 minutes per six-date mapping gives roughly42–55
minutes per 100-date mapping, before accounting for reusable work or a changed
iteration count. This is a planning estimate, not a runtime measurement or a
claim that a full fitted long history will finish overnight. Time a native
long mapping and record memory before granting its subsequent search budget.
Provide progress during backward and forward date sweeps, not only at the end
of a complete long mapping. Independent jobs do not remove the backward/forward
dependence within a forecast or the sequential dependence between shock windows.

Use autonomous Torch controllers and dependency jobs, so laptop sleep does not
stop computation. Checkpoint each completed case and report progress at least
every five minutes. Maintain latest-completed and best-so-far summaries. A
thirty-minute absence of progress is unhealthy and must be diagnosed. Failures
are isolated by case/branch; bounded recovery uses a changed numerical start or
diagnosed repair, not endless identical retries. Collect partial evidence even
when a chain exhausts its budget. No numerical gate is relaxed to meet a deadline.

Use cheap scripted collection and sparse meaningful notifications. Do not
reactivate expensive continuous AI polling or redeem usage credits as part of
this planning request. The formal launch manifest must pin source, objective,
normalization, fiscal, demographic, migration and horizon contracts separately
for A0, A+, B0 and B+. Label each population/entry object as estimated,
empirically normalized, externally fixed or outstanding.

## Existing anchors

- Canonical state and stopped-job record: `CALIBRATION_STATUS.md`.
- Working initial objective: `output/model/e5f_matched_pf_20260909a/initial_calibration_contract/working_contract.json` and `working_weights.csv`.
- Current inspected history/table: `output/model/e5f_matched_pf_20260909a/current_candidate_transition/overnight_20260912/recovered_sequence/`.
- Existing surprise controller: `code/model/tools/run_e5f_successive_surprises_overnight.py`; forecast implementation: `code/model/tools/e5f_successive_surprises.py`.
- Cluster workflow: `docs/workflow/delegation_and_cluster_playbook.md` and `code/cluster/torch.sh`.

These paths are starting references, not approval to reuse stale source or
economic contracts. This document records what to prepare and launch next;
it does not report jobs as submitted.
