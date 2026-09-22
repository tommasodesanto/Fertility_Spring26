# V5 earnings/wealth verified smoke collection

This folder is a mechanical collection of the completed verified smoke outputs from `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/earnings_wealth_direct_period_20260922_v5/output/smoke`. It is numerical smoke evidence only: anchor has two repeated scored evaluations; probe has one scored evaluation and is not an exact repeated run. No refit or paper adoption is implied.

The frozen search initializes from the repeated anchor; the single-repeat smoke probe is outside its selection pool. If the final selected loss exceeds the probe loss, report the better smoke observation explicitly.

Near-bound flags use one percent of the actual raw bound width; they are descriptive and do not establish identification.

Each arm preserves the full 13-row scored target table and 17-row parameter table for every available repetition. The separate `parameters_actual_bounds.csv` applies the frozen V5 structural bounds (including beta upper bound 0.99 and h_P upper bound 2.3) and retains external restrictions for non-structural rows.

- [Collection manifest](manifest.json)
- [Frozen V5 plan](frozen_v5_plan.json)
- [Frozen V5 hash manifest](frozen_v5_hash_manifest.json)
- [Anchor arm](anchor/) — [full target table](anchor/target_fit.csv), [17-row parameters](anchor/parameters.csv), [actual-bound table](anchor/parameters_actual_bounds.csv), [wrapper summary](anchor/wrapper_summary.json), [raw summary](anchor/raw_summary.json), [runtime contract](anchor/runtime_contract.json), and [graphs](anchor/graphs/)
- [Probe arm](probe/) — [full target table](probe/target_fit.csv), [17-row parameters](probe/parameters.csv), [actual-bound table](probe/parameters_actual_bounds.csv), [wrapper summary](probe/wrapper_summary.json), [raw summary](probe/raw_summary.json), [runtime contract](probe/runtime_contract.json), and [graphs](probe/graphs/)

The manifest records source and target SHA-256 hashes, runtime hashes, graph byte/hash verification for all 17 standard graphs per arm, and stationary solve counts from the native raw receipts (six per repetition).

## Standard diagnostic figures

| Figure | Anchor | Inward-beta probe |
|---|---|---|
| fertility by age | [View](anchor/graphs/fertility_by_age.png) | [View](probe/graphs/fertility_by_age.png) |
| fertility policy by age income state | [View](anchor/graphs/fertility_policy_by_age_income_state.png) | [View](probe/graphs/fertility_policy_by_age_income_state.png) |
| housing by age income state | [View](anchor/graphs/housing_by_age_income_state.png) | [View](probe/graphs/housing_by_age_income_state.png) |
| housing market | [View](anchor/graphs/housing_market.png) | [View](probe/graphs/housing_market.png) |
| housing prices | [View](anchor/graphs/housing_prices.png) | [View](probe/graphs/housing_prices.png) |
| income state outcomes | [View](anchor/graphs/income_state_outcomes.png) | [View](probe/graphs/income_state_outcomes.png) |
| liquid wealth by age income state | [View](anchor/graphs/liquid_wealth_by_age_income_state.png) | [View](probe/graphs/liquid_wealth_by_age_income_state.png) |
| market clearing by market | [View](anchor/graphs/market_clearing_by_market.png) | [View](probe/graphs/market_clearing_by_market.png) |
| market clearing residuals | [View](anchor/graphs/market_clearing_residuals.png) | [View](probe/graphs/market_clearing_residuals.png) |
| owner rungs | [View](anchor/graphs/owner_rungs.png) | [View](probe/graphs/owner_rungs.png) |
| ownership by age | [View](anchor/graphs/ownership_by_age.png) | [View](probe/graphs/ownership_by_age.png) |
| ownership by age income state | [View](anchor/graphs/ownership_by_age_income_state.png) | [View](probe/graphs/ownership_by_age_income_state.png) |
| policy childless renter age30 | [View](anchor/graphs/policy_childless_renter_age30.png) | [View](probe/graphs/policy_childless_renter_age30.png) |
| policy childless renter age42 | [View](anchor/graphs/policy_childless_renter_age42.png) | [View](probe/graphs/policy_childless_renter_age42.png) |
| tenure services | [View](anchor/graphs/tenure_services.png) | [View](probe/graphs/tenure_services.png) |
| wealth dist childless renter age30 | [View](anchor/graphs/wealth_dist_childless_renter_age30.png) | [View](probe/graphs/wealth_dist_childless_renter_age30.png) |
| wealth dist childless renter age42 | [View](anchor/graphs/wealth_dist_childless_renter_age42.png) | [View](probe/graphs/wealth_dist_childless_renter_age42.png) |
