# Closed demographic credit benchmark

Verified same-parameter endpoint in `run_v1/`; this is not a verified transition.
Experimental change: remove artificial purchase, collateral and unsecured credit limits, retaining feasible continuation and nonnegative net estates at every possible death date. The approximation uses feasible grid nodes and the inherited value cutoff. All 31 parameter estimates, earnings, conditional entrant wealth, preferences, grid, weights and targets remain fixed. In particular, child utility is not recalibrated.

A closed stationary population requires births to replace entry. Price solves that renewal condition; population then scales to the unchanged absolute housing supply curve. This is why long-run completed fertility returns to 2.100: population and prices, rather than preference renormalization, adjust.

| Outcome | Reference | Closed credit endpoint |
|---|---:|---:|
| House price | 0.826 | 0.852 |
| Population | 1.000 | 1.035 |
| Completed fertility | 2.100 | 2.100 |
| Ownership ages 30–55 | 0.617 | 0.635 |
| Wealth / annual earnings | 6.025 | 5.640 |
| Early fertility | 0.528 | 0.553 |
| Common weighted loss | 42.282 | 57.193 |

Price rises 3.133% and population 3.524%. The benchmark improves ownership and early fertility but leaves early fertility well below its 0.810 target and worsens wealth fit. It is a counterfactual, not a recalibration or a better-fitting candidate.

All 14 target rows, targets, gaps, weights and contributions: `run_v1/target_comparison.csv` and `run_v1/selected/target_fit.csv`. All 31 parameters including the ten fitted reference parameters, their original bounds and bound flags: `run_v1/selected/parameters.csv`. The fertility normalization row is a restriction, not a weighted loss row.

Twelve fixed-price solves including control and exact repeat completed within the bounded plan. Renewal residual is -1.810e-7; absolute housing demand equals supply. All 14 target and 31 parameter repeat rows are identical. Estate funding passes without negative estates. The 17 standard figures are in `run_v1/selected/standard_diagnostics/`, with an inspected contact sheet at `run_v1/diagnostic_contact_sheet.png`. Control and repeat receipts remain in the run folder.

The selected checkpoint contains a UNIT-MASS stationary distribution. A transition terminal-distribution comparison must multiply it by the recorded population scale 1.0352351249517946. Value functions require no population scaling. Never overwrite the forward-carried transition distribution with this endpoint.

Reproduction driver: `code/model/tools/run_e5f_credit_closed_endpoint.py`; source and runtime pins in `plan_v1.json`. Completed outputs remain immutable. Production native replay and full transition certification are separate work indexed by `../credit_transition/preparation/PRODUCTION_INTEGRATION.md`.
