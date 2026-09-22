# V5 terminal search evidence

The search stopped on proposal 006's native timeout before final repeated verification. Proposal 005 is the best of seven scored proposals out of eight attempts: loss 1087.2816148857435, evaluated once. The earlier single-repeat smoke test 1040.4468628084285 remains better and was outside the frozen search selection pool. No adopted or finally verified calibration is claimed.

- [Reviewed economic and numerical assessment](../terminal_review.md)
- [Full 13-row target table](target_fit/target_fit_case005.csv)
- [All 17 parameters with actual bounds and near-bound flags](actual_bounds/parameters_actual_bounds_case005.csv)
- [Raw parameter table](parameters/parameters_raw_case005.csv) and [scored parameter table](parameters/parameters_scored_case005.csv)
- [All 8 attempted proposals and honest solve counts](search/attempted_proposals.csv)
- [Original scored-case stream](search/cases.jsonl)
- [Source, target and runtime fingerprints](receipts/runtime_contract.json), [wrapper receipt](receipts/wrapper_summary.json), and [score](receipts/score.json)
- [Native timeout](receipts/case006_failure.json), [seventh-solve heartbeat](receipts/case006_raw_heartbeat.json), [controller failure](receipts/search_failure.json), and [terminal supervisor receipt](receipts/execution.json)
- [Immutable V5 plan](../staging/frozen_v5_plan.json), [source inventory](../staging/frozen_v5_hash_manifest.json), and [collection hash manifest](manifest.json)

Search counts: 49 started,48 completed, 1 incomplete stationary solve. V5 including smoke: 67/66/1; all prior attempts included: 102/98/4. Counts distinguish completed stationary solves from complete scored objectives. There were seven scored search objectives, ten scored V5 objectives including smoke, and twelve including the two historical V4 anchors. No model execution or checkpoint replay occurred during collection.

Near-bound flags use 1% of the actual raw bound width. Annual beta's upper bound is 0.99; h_P's is 2.3. Generic scorer bounds in raw tables remain preserved; use the annotated table for the actual experiment.

## Original standard figures

All 17 figure bytes were independently checked against the wrapper hashes. The stable set is unchanged; 45-state legends and the long wealth axis limit readability.

- [fertility by age](standard_diagnostics/fertility_by_age.png)
- [fertility policy by age income state](standard_diagnostics/fertility_policy_by_age_income_state.png)
- [housing by age income state](standard_diagnostics/housing_by_age_income_state.png)
- [housing market](standard_diagnostics/housing_market.png)
- [housing prices](standard_diagnostics/housing_prices.png)
- [income state outcomes](standard_diagnostics/income_state_outcomes.png)
- [liquid wealth by age income state](standard_diagnostics/liquid_wealth_by_age_income_state.png)
- [market clearing by market](standard_diagnostics/market_clearing_by_market.png)
- [market clearing residuals](standard_diagnostics/market_clearing_residuals.png)
- [owner rungs](standard_diagnostics/owner_rungs.png)
- [ownership by age](standard_diagnostics/ownership_by_age.png)
- [ownership by age income state](standard_diagnostics/ownership_by_age_income_state.png)
- [policy childless renter age30](standard_diagnostics/policy_childless_renter_age30.png)
- [policy childless renter age42](standard_diagnostics/policy_childless_renter_age42.png)
- [tenure services](standard_diagnostics/tenure_services.png)
- [wealth dist childless renter age30](standard_diagnostics/wealth_dist_childless_renter_age30.png)
- [wealth dist childless renter age42](standard_diagnostics/wealth_dist_childless_renter_age42.png)
