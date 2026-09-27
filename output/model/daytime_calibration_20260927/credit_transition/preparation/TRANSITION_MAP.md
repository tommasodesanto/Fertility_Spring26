# Read-only transition mapping — September 27

Scope: exploratory permanent removal of artificial credit limits, fixed preferences, initialized at verified overnight de_0093. No model run, core edit, or economic closure adopted in this mapping.

## Reusable implementation

The frozen September 27 source under `tmp/e5f_overnight_local_20260927/portable/calibration_code_integration_20260927_v2/source/code/model/` includes `tools/run_e5f_perfect_foresight_transition.py`. Prefer importing its functions through the authenticated runtime already used by `run_e5f_credit_benchmark.py`, not executing its historical CLI: CLI defaults reconstruct an August 2023 economy and solve an unrelated preference endpoint.

Useful functions in that module:
- `solve_date_policy` (line488): invokes the authenticated Bellman with `continuation_V`.
- `backward_value_path` (508): backwards calendar recursion, takes fixed psi path, terminal value, price/rent, pension and payroll-tax paths.
- `evaluate_path_at_prices` (561): forward level-valued household distribution, child birth queue, housing demand and mass audits. Takes caller-supplied supply and birth-to-entry conversion, so no need to adopt old defaults.
- `rents_from_asset_prices` (297): exact stationary-compatible dated user cost, r_t=(R+delta+tau_H)p_t-p_{t+1}; rejects nonpositive rents.
- `terminal_convergence_diagnostics` (346): compares actual carried terminal population/distribution/entry queue with endpoint, rather than replacing the state.

For a joint housing/PAYGO path use the generic, caller-owned `e5f_social_security_root.solve_social_security_path`, with current baseline payroll tax held fixed and pensions endogenous to actual dated workers/retirees. It supports evaluation cap, absolute deadline, callback, fresh replay and `initial_jacobian`/returned `final_jacobian`. Do not copy the .179 tax hardcoded in older drivers: the September27 runtime derives its own baseline tax. Property-tax rebates are a separate closure; do not introduce them merely because an old driver was named rebated.

`e5f_original_queue_experiment.queue_path` demonstrates wiring the existing forward law, but hardcodes old payroll tax and birth conversion and hence is not drop-in. `run_e5f_final_rebated_history.solve_forecast` demonstrates dated audit/heartbeat/root wiring; it also carries old rebated fiscal and person/queue assumptions, so use as a pattern only.

## Credit adapter work before transition

`e5f_solvency_credit_benchmark.install` currently explicitly REJECTS non-None `continuation_V` (line95). It is stationary-only by design. A fresh version must allow the dated boundary and regression-test it, leaving completed run_v2 immutable.

The underlying frozen Bellman already chooses `next_values = V if continuation_V is None else continuation_V` (solver.py:2578). Thus the strict positive-probability infeasibility propagation patch can naturally inspect tomorrow's value tensor. This does not alone certify the whole transition: test stationary constant paths reproduce existing policies, probabilities, 14 targets and fiscal accounts; test next-date bankruptcy remains infeasible for every positive income/child transition. Preserve the explicit V>-1e9 classification and grid-node approximation limitation.

Timing needs an explicit check: the new death-solvency mask uses `b' + (1-selling_cost)*p*H`, where the stationary Bellman builds bequest house value at its passed price. On a nonconstant path, verify which dated house price values assets at death and reconcile that with existing owner/renter no-arbitrage, continuation assets, and the independent estate ledger. Do not simply drop the continuation guard and assert success.

Initial states must use de_0093 pre-fertility distribution, grid, parameters and raw/adjusted entry pipeline; not post-choice distribution or the older reconstructed2023 distribution. Build a constant-policy baseline path first to prove the initial queue timing and mass law reproduce its steady state.

## Demographic/fiscal/supply closure decisions

A conditional unit-mass stationary cross section is not closed renewal. The new endpoint must hold psi fixed and allow population scale and price to adjust jointly to renewal, while balancing PAYGO. Holding population fixed and re-normalizing entry each date would undo the requested experiment. A fixed-preference positive stationary endpoint may fail to exist; never tune psi or migration silently to obtain one.

The existing four-vintage queue tracks both raw and top-code-adjusted births and advances unnormalized household mass. It historically uses births/2.1 and a sixteen-year maturation delay. Whether this is the chosen September27 diagnostic demographic closure must be recorded explicitly, especially with model entrants beginning at18. The more elaborate person-cohort code is an additional economic model, not an automatic replacement. Source helpers permit an explicit chosen conversion and queue, so the lead can use the intended current contract without adopting historical constants.

Use the same absolute housing supply schedule as de_0093 (elasticity .63, inherited H0/user-cost normalization), not a newly per-capita-rescaled supply after population changes. Carry current entry wealth/income distribution, fixed interest, no migration if a closed experiment is intended, and current no-rebate property-tax convention. At each date record positive/negative estates, entry wealth costs, residual estate sink or shortfall; no donor insurance/default or transfers silently added. Preserve gross versus net bequest valuation convention.

## Speed/cost: test before long launch

The native path mapping solves each date twice: H backward Bellmans plus H forward policy reconstructions (`backward_value_path`, then `solve_date_policy` in the forward loop). Thus K price/fiscal root evaluations cost roughly 2HK Bellmans plus terminal roots, audits and replay. H=8,K=10 is160 Bellmans; H=16,K=20 is640. These are full lifecycle Bellmans, not single-age updates.

Latest relevant measured evidence: credit run_v2 natural stationary solve70.148s,4 price evaluations (`natural/case/stationary_solves.json`), roughly17.5s per equilibrium price evaluation including stationary distribution/overheads. This is NOT a measured dated-Bellman timing. If a dated call cost a similar amount, H8/K10 would be about47minutes before endpoint/export; profile H2 first instead of promising the transition is as fast as one SS.

Low-risk existing acceleration: `e5f_exact_policy_cache.policy_cache(pf,max_bytes=...)` hashes complete call state and can reuse identical backward policies in forward reconstruction. It has bounded memory and hit/miss/eviction/serialization-bypass counters. Prove identical arrays with cache off/on; full policy storage can be large, and eviction order can wipe the hoped-for hits. Reusing the backwards bundles directly is a possible cleaner change, but needs explicit memory accounting and separate regression review.

Existing root improvement: supply a measured small directional/structured initial Jacobian and retain subsequent root Jacobian warm starts. Historical `docs/model/e5f_sequence_space_prototype.md` reports a different ten-date problem around3min/mapping and a7-mapping block-Toeplitz initialization reducing8-update residual1.346 to.0225. These are evidence of a useful method, not current performance promises or a ready-made derivative for this new credit regime. Avoid expensive dense 2n finite-difference construction before profiling. The local calibration Jacobian is a moments-vs-parameters derivative, not the dated price/PAYGO Jacobian needed here.

Recommended bounded first pass: exact constant-path H2 control; H2 credit mapping and solvency/estate timing checks; cache on/off equality plus phase times; then a versioned H8 exploratory path with explicit root/evaluation/time caps and checkpoint per mapping/date. Label finite-horizon terminal gaps honestly. Extend horizon only after the endpoint and short mapping are verified.

No new numerical result is claimed here. Lead owns endpoint closure, timing resolution and launch budget.

## Author-requested September14 reference reconciliation (supersedes generic closure discussion)

`tmp/paper_baseline_sep14/PAPER_BASELINE.md` is authoritative for the preserved presentation calculation: balanced PAYGO, tax .179, original household birth queue with2.1 replacement, no immigration,1%annual property tax with equal rebates. Its actual presentation path was four announced shocks and unconverged, not a certified single permanent shock. These are historical source facts; reconcile explicitly with September27 changes rather than importing every historical numeric constant.

Exact queue initialization is documented/implemented in `tmp/paper_baseline_sep14/output/model/paper_baseline_sep14/native/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/afternoon_original_queue_20260913a/source/e5f_original_queue_experiment.py:44`: preserve exact stationary pre-fertility tensor; initialize four adjusted waiting slots with actual age-zero mass; initialize four raw slots with raw `total_births_kfe/2.1`. No person-state reset, age reweighting, immigration or split-birth queue. IMPORTANT correction to the earlier generic description: four waiting slots enter the forward state at t+5, a20year birth-to-entry effect, NOT16years. `run_e5f_post2023_no_policy_continuations.py` metadata explicitly reports waiting_slots+1 dates.

Within-period historical timing: `tmp/paper_baseline_sep14/docs/model/timing_repair_advisor_recap_20260723.md:24` records beginning assets/tenure -> choices and housing transactions -> posttransaction liquid wealth -> consumption/saving b' -> survival/death. The estate statistic is post-saving b'+pH of newly chosen tenure; living wealth remains beginning b+pH of inherited tenure. Frozen Bellman and statistic use the passed/current price. This supports preserving current-price death valuation unless an explicit nonstationary decision establishes otherwise. The recap itself has no t-versus-t+1 price indexing, so it does not mathematically settle all nonstationary pricing conventions.

Do not import `tmp/paper_baseline_sep14/docs/model/simplified_olg_transition_max_theory_audit_20260830.md`'s stepped-up estate valued at next-period price into this quantitative solver: it analyzes the distinct two-generation theory. The historical quantitative PF code already combines the current-price death table with next-date continuation values and its exact rental identity. Report that as the inherited convention, and test it consistently rather than change it while calling the exercise frictionless credit only.

The native `run_e5f_post2023_no_policy_continuations.solve_closed_stationary_endpoint` already implements fixed-psi closed renewal: root house price on adjusted births/entry=2.1, then scale population as absolute supply divided by per-household housing demand. This is the correct reusable endpoint architecture for the author's intended inherited closure; no need to construct a person/head replacement. Pension balance and September27 credit/utility runtime still require current adapters.

Implementation remains paused per lead's historical-reconciliation request. No experiment, adapter edit or launch has been made by this mapping task.
