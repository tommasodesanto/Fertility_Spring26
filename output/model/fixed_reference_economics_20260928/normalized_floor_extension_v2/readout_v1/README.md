# Normalized floor extension v2 readout

Collected the four completed diagnostic points from Torch job array `19002367`, preserving each remote `completed.json`, launcher start/terminal receipts, native ROOT and REPEAT target fit, parameter table, closure JSON, and all 17 standard diagnostic PNGs. No caches, arrays, checkpoints, or model reruns were included.

`summary.csv` retains every column from the native 14-row target fit and 31-row parameter table for each point. It marks the nine incumbent search coordinates as `fixed_at_incumbent`, `h_P` as `diagnostic_varied`, H0 as `derived`, and remaining retained parameters as `fixed_in_diagnostic` (their native status columns retain their original economic roles). The full row-level loss contributions, weights, estimates, bounds, and near-bound fields are carried through from the native CSVs.

The points are h_P = 2.3, 2.4, 2.5, and 2.6, with recorded losses 29.9519, 42.1673, 79.8038, and 149.1804. This readout reports run outputs only and makes no economic interpretation.

Every collected file's SHA-256 matches the remote collector manifest. For every point, ROOT and REPEAT target CSV, parameter CSV, closure JSON, and all 17 PNGs are byte-identical. All four launcher receipts record exit 0 and matching Slurm job IDs between start and terminal receipts.

The 17 standard graph files are: ``fertility_by_age.png`, `fertility_policy_by_age_income_state.png`, `housing_by_age_income_state.png`, `housing_market.png`, `housing_prices.png`, `income_state_outcomes.png`, `liquid_wealth_by_age_income_state.png`, `market_clearing_by_market.png`, `market_clearing_residuals.png`, `owner_rungs.png`, `ownership_by_age.png`, `ownership_by_age_income_state.png`, `policy_childless_renter_age30.png`, `policy_childless_renter_age42.png`, `tenure_services.png`, `wealth_dist_childless_renter_age30.png`, `wealth_dist_childless_renter_age42.png``. The PNGs are under each `pointN/ROOT/standard_diagnostics/` and `pointN/REPEAT/standard_diagnostics/` folder for immediate visual inspection. See `verification.json` for hashes, row counts, job IDs, and check results.
