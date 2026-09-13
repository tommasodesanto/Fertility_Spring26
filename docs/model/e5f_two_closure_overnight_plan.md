# Two demographic closures: overnight rerun plan

**Author decision, September 12–13, 2026. Planning record only; no launch from this request.**
Run both demographic specifications. Prefer the survivor-maturation version if
it produces verified, usable results; retain the current specification as the
fallback, with its missing treatment of orphan care stated explicitly. Neither
branch is selected merely because its job finishes or the other branch fails.

## The two branches

| Branch | Demographic treatment | Presentation role |
|---|---|---|
| A: retained demographic treatment | Keep today's historical head-age conditioning, 2023 person anchor and subsequent person-cohort/headship law. The household operator removes dependents when the parental household exits, while the separate population calculation preserves its own surviving child cohorts. No explicit consumption or housing expenditure is assigned to these orphans. | Fallback, with this omission disclosed; do not describe it as an implemented or financed care system. |
| B: entry from surviving maturation | Dependents die when their parental household exits. Apply the existing stochastic maturation process only to surviving dependents. Domestic household entry comes only from that surviving maturation flow, with an explicit conversion/formation rule and separately counted migration. | Preferred candidate, conditional on numerical and empirical checks. Joint family mortality is an assumption, not a demonstrated empirical improvement. |

Keep source snapshots, output folders, target receipts and checkpoint chains
separate. Never splice A's initial state, demographic anchor or terminal
equilibrium into B's history. In B, the separate birth queue and person/head
rescaling cannot silently recreate entrants or overwrite the new closure.

## Shared economic specification

- Use the current sequential-choice model, revised utility, stochastic child
  maturation and nine-parameter structural search. Do not reopen nesting or
  utility experiments tonight.
- Estimate annual beta with its existing cap of 0.99; do not fix it at the cap.
  Preserve the remaining approved bounds and external restrictions.
- **Rebate all property-tax revenue equally in both baselines.** Re-solve the
  initial equilibrium, every forecast and the terminal equilibrium with this
  same fiscal rule. PAYGO pensions balance separately at each date.
- Refit the approximate pre-2007 stationary economy using the existing twelve
  scored empirical moments, their definitions and weights. Keep the first-birth
  room response at 0.7202462623815278. No target dropping or reweighting.
- Fit the four fertility windows from 2007 through 2023 using successive
  unexpected preference changes. At each change, households expect that level
  to persist; carry the realized household state forward. Hold preferences
  constant after 2023. These are fitted deterministic surprises, not an
  estimated stochastic preference process.
- The complete 2023 Data/Model table remains an untargeted validation exercise.
  Preserve its empirical-vintage labels and all thirteen moment families.

## Launch sequence and parallel work

### 1. Two exact-loop smoke jobs, before large searches

For A, verify the equal-rebate initial, terminal and dated fiscal calculations.
For B, first exercise the demographic operator at saved policies, then connect
it to the same fiscal and equilibrium solver. No additional household states
are required for an aggregate stochastic-renewal experiment.

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

The old initial 2.1 normalization is not proof of replacement under B. Check
actual reproduction: stationary household entry must equal household exits.
Preserve A's current normalization for its controlled rerun. If B needs a
different normalization or formation restriction, record that explicitly before
its search; do not force both conditions by undocumented rescaling. Its twelve
scored empirical rows remain unchanged, but its normalization contract must be
separately named. This is the remaining specification item to settle at launch,
not a reason to leave A idle.

Both exact loops must write checkpoints, full observations, failed-case
receipts and the standard diagnostic packet before a search is admitted.

### 2. Full initial calibration in both branches, in parallel

Use all nine estimated coordinates and the complete objective. Proposed search
ceiling per branch: eighteen plus/minus coordinate probes, then up to two
eighteen-candidate waves of joint proposals informed by the complete objective
and local response matrix. These joint waves may move all nine parameters.
Keep the best valid point and perform two independent exact repetitions.
No branch wins from a favorable subset of moments.

Use up to eighteen concurrent calibration workers in total, subject to actual
Torch memory/account limits; each worker has one numerical thread. Start a
historical pipeline from the first verified rebated initial candidate in each
branch while the structural search continues. Freeze its source and parameters.
A better initial calibration can start a separate history; it cannot silently
replace the starting distribution of an existing fitted path.

### 3. Historical shock fitting and horizon checks

Give each branch up to two independent historical searches with different
numerical starts. Within each history, windows must run sequentially because
the next window inherits the preceding realized state. Use bracketed searches
over preference levels, saved price/pension/transfer guesses and exact replay.
Reject a failed trial and continue with the next admissible proposal; do not
terminate unrelated chains or admit a failed inherited state.

Six-date forecasts are initialization diagnostics. Compare successively longer
forecasts (initially twelve and twenty-four four-year dates), preserving the
same inherited state and economic specification. Refit shock levels if changing
the horizon materially changes the historical fertility fit. Extend further
only within the measured run budget. An endpoint-distance failure or material
change with the horizon stays visible; a longer attempted run is not a
certificate. Retain the existing fertility tolerance and numerical gates.

### 4. Policy computations in both branches

The priority comparison is the baseline 1% annual property tax with equal
rebate versus a 2% annual tax with equal rebate. First compute consistent
stationary endpoints in parallel with historical work. Label these long-run
comparisons; they are not effects from the inherited 2023 economy.

As soon as a branch has an admissible fitted 2023 state and baseline
continuation, run its rebated-tax policy transition from that same state.
Compare births, population under that branch's definition, housing services,
ownership, prices/rents, consumption and both fiscal budgets. Baseline and
policy require matching horizons and demographic closure.

Supply +20% and dependent-child LTV95% are the next independent policy jobs if
the main historical/tax work passes and time remains. Policy failure must not
erase completed calibration, history or stationary evidence. The present
six-date, horizon-unverified path cannot be silently promoted to a production
policy benchmark.

### 5. Automatic collection and morning decision

Write one review packet per branch and a side-by-side comparison containing:

- every initial target, model moment, gap, weight and loss contribution;
- every estimated parameter, bound and bound proximity, plus external inputs;
- all four fertility windows, fitted preferences, model/data path and saved
  continuation, with clear forecast-horizon status;
- all thirteen 2023 Data/Model rows with source vintages;
- the stable policy-function, lifecycle, intergenerational-allocation, price,
  quantity, market, population and fiscal diagnostic graphs;
- policy differences from a matched baseline, separating stationary evidence
  from transition evidence, and a concise list of failed or outstanding gates.

Prefer B only if its accounting, reproduction, household choices, fiscal and
market checks pass and its fit and horizon evidence support the stated claims.
Otherwise use A only to the extent its own checks permit, explicitly stating
the orphan-care omission. If neither passes a transition/horizon gate, show
provisional results as provisional. Do not imply that either specification is
globally unique or guaranteed to exist.

## Compute budget and unattended operation

The proposed envelope is twelve hours from launch, with the final hour reserved
for reproduction, collection and figures. Initial search should use at most
three hours; it must not delay the first viable historical pipeline. Historical
chains have explicit trial and per-forecast limits and a shared deadline;
reserve the final three hours for admitted policy runs and verification. These
are caps, not forecasts of successful completion.

Search size ceiling: 54 objective evaluations plus two smoke repetitions and
two selected-point repetitions per branch, or 116 single-repetition evaluations
across both. With at most eight stationary solves per such evaluation, this
means at most 928 stationary solves. At the previously observed 2–5 minutes per
stationary solve, that is roughly 31–77 CPU-hours, or 1.7–4.3 hours with eighteen
fully utilized workers, before startup, uneven work, failures and queueing.
Replace this rough estimate with the exact-loop measured cost before submission;
reduce round counts if the three-hour stage budget requires it.

Historical fitting has up to four chains, four windows and six preference
trials per window: at most 96 forecast attempts before additional horizon
checks, with time caps usually binding much earlier. A six-date converged
forecast previously took 6–45 minutes at fixed preferences; longer forecasts
and repeated shock roots can dominate total time. Parallelism helps independent
chains and policy branches, not the dependence between historical windows.
Before launch, record horizon-specific solve counts, memory and wall-time
estimates from the native smoke; no uncosted 24-date or longer search.

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
normalization, fiscal, demographic and horizon contracts separately for A and B.

## Existing anchors

- Canonical state and stopped-job record: `CALIBRATION_STATUS.md`.
- Working initial objective: `output/model/e5f_matched_pf_20260909a/initial_calibration_contract/working_contract.json` and `working_weights.csv`.
- Current inspected history/table: `output/model/e5f_matched_pf_20260909a/current_candidate_transition/overnight_20260912/recovered_sequence/`.
- Existing surprise controller: `code/model/tools/run_e5f_successive_surprises_overnight.py`; forecast implementation: `code/model/tools/e5f_successive_surprises.py`.
- Cluster workflow: `docs/workflow/delegation_and_cluster_playbook.md` and `code/cluster/torch.sh`.

These paths are starting references, not approval to reuse stale source or
economic contracts. This document records what to prepare and launch next;
it does not report jobs as submitted.
