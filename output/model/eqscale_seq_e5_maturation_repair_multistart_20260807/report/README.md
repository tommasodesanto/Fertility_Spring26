# E1 collector

This collector scans `production/chain_*/summary.json`.
It retains only strict, exactly repeated tight winners.
Exact equality is required for loss and target moments.
The lowest eligible tight-repeat loss is the E1 winner.
`results.json` exposes `winners.E1` in the M5 winner shape.
`results.json` also preserves the selected chain's complete calibration metadata.
`results.json` preserves the exact-repeat check and source summary path.
`target_fit_full.csv` contains every row in the selected target system.
`parameter_table_full.csv` contains every estimated parameter, its search bounds, and a two-percent near-bound flag.
`chain_summary.csv` records eligibility for every discovered chain.
No model solve is performed by this collector.
