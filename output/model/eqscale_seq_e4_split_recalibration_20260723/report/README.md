# E1 collector

This collector scans `production/chain_*/summary.json`.
It retains only strict, exactly repeated tight winners.
Exact equality is required for loss and target moments.
The lowest eligible tight-repeat loss is the E1 winner.
`results.json` exposes `winners.E1` in the M5 winner shape.
`target_fit_full.csv` contains all 15 unchanged target rows.
`chain_summary.csv` records eligibility for every discovered chain.
No model solve is performed by this collector.
