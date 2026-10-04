# Two unanticipated fertility-preference shocks

**Status:** preparation only, October 4, 2026. No native two-shock estimation has launched. The [one-shock experiment is retained separately](../retained_one_shock_v1/README.md).

The author requested two fitted anchors: retain the final 2020–2023 fertility target and add a midpoint target. The working implementation uses the existing **2012–2015** midpoint block and shocks in **2007 and 2015**. In each period households believe the current preference level is permanent; the 2015 change is unanticipated in earlier decisions.

| Birth window | Decision year | Target | Weight | Role |
|---|---:|---:|---:|---|
| 2008–2011 | 2007 | 1.974875 | 0 | Validation |
| 2012–2015 | 2011 | 1.861000 | 1 | Fit first shock |
| 2016–2019 | 2015 | 1.755375 | 0 | Validation |
| 2020–2023 | 2019 | 1.645750 | 1 | Fit second shock |

The two free parameters are the permanent preference levels after each surprise. Both retain baseline **0.17892072066041628** times **[0.01, 2]**, or absolute bounds **[0.001789207206604163, 0.35784144132083257]**. Neither is estimated yet. Two informative fitted moments are required, and each scalar fit must demonstrate a nonzero local response and freshly reproduce its accepted equilibrium.

First fit the 2007 surprise to local fertility index 1. Replay only its accepted first two periods, using that forecast's own 2015 price and value function as the boundary. Save the complete inherited 2015 household distribution and both adjusted/raw entry queues; check exact state transfer and strict prefix replay. Then reveal the second surprise and fit its local index 1 to the final target. The 2015 boundary price and pension are numerical initial guesses for the second equilibrium, not fixed prices. Earlier choices must not use the second surprise's forecast.

**Changes relative to the retained one-shock experiment:** the second preference surprise and promotion of the existing midpoint validation row to a fitted row are author-requested experimental changes. The first shock is re-estimated. Target definitions, sample, four-year clock and target values remain unchanged. Earnings, wealth/income distributions, post-interest transaction timing, soft financing, net estates, housing supply, entry rules, both birth queues, fixed payroll tax/endogenous pension, and all baseline parameters are unchanged. Keep saved one-birth Estate-A case `20261003T212605039706Z_a739edc3`; do not adopt concurrent recalibration or birth-cap extensions.

The numerical design reuses the scalar fitter and native mappings in a separate controller. Each stage has distinct state-bound numerical guesses and a fresh five-map, 12-date Jacobian. The second seed is measured at the actual inherited 2015 state; stationarity is tested at its stationary endpoint, not imposed on the inherited distribution. Both 24/32-date paths must pass the original root, accounting, replay and historical-horizon checks. Terminal failures remain explicit. This diagnostic does not certify 104/128-period production paths or close the provisional estate contract.

Before a long fit: authenticate the isolated source package; review surprise timing, target indices and both queue transfers; pass focused tests and an exact-loop native smoke; verify checkpoint, latest/best summaries and stable diagnostic plots. Proposed local allocation is one core, 24 GiB, 21,480 internal seconds inside six hours, 20,000 actual native calls shared across both stages, 12 fit evaluations per stage, 12 path evaluations, and 48 stationary endpoint iterations within 1,800 seconds. Cost and the exact launch receipt must be recorded before launch. No old run is resumed or extended.
