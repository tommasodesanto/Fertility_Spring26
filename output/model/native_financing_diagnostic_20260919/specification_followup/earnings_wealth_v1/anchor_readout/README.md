# Earnings-wealth V4 anchor readout

This is a bounded collection of the V4 smoke anchor, status `verified_scored_candidate`. It records the starting parameter vector and two exact scored repetitions at loss `1160.3761384370134`; it is not a refit, convergence claim, or adopted calibration.

The complete 13-target raw fit is [`target_fit.csv`](target_fit.csv), and the complete 17-parameter raw table is [`parameters.csv`](parameters.csv). The second exact repetition is retained under [`repetition_02/`](repetition_02/). [`parameters_actual_bounds.csv`](parameters_actual_bounds.csv) annotates all 17 parameters while applying the plan's actual bounds only to the nine structural coordinates; the generic raw bounds remain unchanged in `parameters.csv`. Near-bound means within 1% of the actual structural-bound range. The actual upper bounds are beta `0.99` and `h_P` `2.3`.

The native wrapper summary, raw native summary, score, run contract, runtime contract, and source/objective fingerprints are retained in this directory. The native run used 12 stationary solves across two repetitions (six each); both repetitions have 13 target rows and 17 parameter rows, and their losses are exactly equal.

The 17 original standard diagnostic graphs are in [`graphs/`](graphs/), copied from `summary.original_graphs` with byte hashes verified. No supplemental graph was added. The age-30 policy display has a crowded legend and a wealth axis extending to 3000; this is a display limitation only, with no economic conclusion drawn here.

The full target rows agree exactly across repetitions apart from the checkpoint-file hash; raw parameter tables are byte-identical. Native assertions verified exact price, value function, distribution, moments and normalized fertility utility across repetitions.

## Standard figures

- [fertility by age](graphs/fertility_by_age.png)
- [fertility policy by age income state](graphs/fertility_policy_by_age_income_state.png)
- [housing by age income state](graphs/housing_by_age_income_state.png)
- [housing market](graphs/housing_market.png)
- [housing prices](graphs/housing_prices.png)
- [income state outcomes](graphs/income_state_outcomes.png)
- [liquid wealth by age income state](graphs/liquid_wealth_by_age_income_state.png)
- [market clearing by market](graphs/market_clearing_by_market.png)
- [market clearing residuals](graphs/market_clearing_residuals.png)
- [owner rungs](graphs/owner_rungs.png)
- [ownership by age](graphs/ownership_by_age.png)
- [ownership by age income state](graphs/ownership_by_age_income_state.png)
- [policy childless renter age30](graphs/policy_childless_renter_age30.png)
- [policy childless renter age42](graphs/policy_childless_renter_age42.png)
- [tenure services](graphs/tenure_services.png)
- [wealth dist childless renter age30](graphs/wealth_dist_childless_renter_age30.png)
- [wealth dist childless renter age42](graphs/wealth_dist_childless_renter_age42.png)
