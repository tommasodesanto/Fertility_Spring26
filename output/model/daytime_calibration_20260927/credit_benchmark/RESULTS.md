# Verified first-pass borrowing comparison

All three cases completed. Baseline replay reproduces all31 parameter estimates
exactly and all14 moments within2.34e-6 of the authenticated overnight solution.
The benchmark and independent repeat reproduce every target and parameter CSV
byte exactly. All economic parameters are fixed, including child benefit.
This is an experimental credit change, not a new calibration.

| Outcome | Reference | Solvency-credit benchmark | Data target |
|---|---:|---:|---:|
| Completed fertility | 2.1001 | 2.1283 | 2.1000 |
| Children by age25 (capped3) | 0.5275 | 0.5642 | 0.8095 |
| Ownership ages30–55 | 0.6174 | 0.6465 | 0.6763 |
| Wealth / annual earnings | 6.0247 | 5.6474 | 6.9266 |
| First-birth housing response | 1.6099 | 1.7132 | 1.4650 |
| Common weighted fit score | 42.282 | 66.439 | — |

Credit relaxation advances fertility and ownership somewhat but does not close
the early-fertility gap; wealth and first-birth housing fit worsen. These are
same-parameter equilibrium comparisons, not identified causal estimates for data.
All14 target values, gaps, weights and contributions appear in
`run_v2/target_comparison.csv`. All31 parameter estimates, original bounds and
bound flags appear in `run_v2/natural/case/parameters.csv`; the10 previously fitted
parameters are now fixed. Full receipts are retained beside each case.

The benchmark steady state takes70.15seconds (repeat70.00); market residual
5.66e-6. Household budgets, full distribution operator, choice probabilities,
occupied value monotonicity and pension balance pass. No occupied transaction
or saving lies outside the grid; no saving mass at either endpoint. Net-negative
estate exposure, death mass and liabilities are exactly zero. Positive net
estates0.102184 fund entrant wealth0.031700, leaving sink0.070485 per model period.
All17 standard plots were inspected using the retained contact sheet. Late-life
ownership still approaches one: numerical success does not remove that economic
miss. Standard policy panels span the full grid and cannot substitute for zoomed
inspection of the populated wealth region.

Important scope: the grid boundary is conservative and uses the native feasibility
value cutoff; this is not a continuum frictionless-credit proof. Solvency at death
is an explicit experimental restriction and is not guaranteed to relax every
unoccupied baseline state. Normalized entrants are inherited, not adjusted to
new births; birth/entry replacement differs by1.33%. Thus this is a conditional
stationary population comparison, not a fully closed demographic counterfactual.
The original preference, estate valuation and funding conventions are unchanged.

run_v1's wrapper error and original sources remain preserved; see README.md.
No further benchmark jobs are running. Cluster search and Jacobian are separate.
