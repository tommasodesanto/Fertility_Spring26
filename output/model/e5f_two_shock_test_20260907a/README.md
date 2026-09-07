# Two contemporaneous shocks: computational test

Completed 7 September 2026. Experimental source: branch
`codex/two-contemporaneous-shocks`, numerical implementation at `910a9464`.
No production changes; no new calibration or market-clearing claim.

**Result:** the operator and full lifecycle test pass. At unchanged parameters
and inherited prices, births change by 0.009572% and ownership by
-0.681636 percentage points.

| Fixed-price lifecycle quantity | Old specification | Two shocks |
|---|---:|---:|
| Births per four-year period | 0.077453471072 | 0.077460885274 |
| Ownership share | 0.596200166097 | 0.589383802880 |
| Mean parity | 1.254643363737 | 1.254763463215 |
| Total population mass | 1.000000000000 | 1.000000000000 |
| Relative housing-market residual | 0.269951036 | 0.266541755 |
| Solve seconds | 16.463 | 111.025 |

These ownership figures cover all adult model households; they are not the
age-restricted calibration ownership moment. Births are period flows, not
completed fertility. Both cases use freshly computed lifecycle continuation;
the historical checkpoint supplies parameters and prices, not a stationary
market-clearing price for these fresh populations. The roughly 27% housing
residual is reported explicitly; this experiment is not an equilibrium.

The two differences are independent centered logistic shocks with housing scale
0.005, first-birth scale 2.1681730392479377 and later-birth scale
1.7364706586958831. No ordering constraint and no GEV lambda. The plan menu
retains committed tenure across conception outcomes and deterministic product
choice within tenure, so the old/new comparison includes those menu differences.
The poor overnight nested-GEV calibration therefore does not establish a poor
fit for this different two-shock specification.

Verification:

- Sixteen tests pass (`tests.log`), including all feasibility masks, extreme
  scale ratios, independent numerical integration and value derivatives.
- All ten default-off arrays reproduce the pristine parent exactly in the same
  local runtime (`final_comparison/baseline_reference.json`). The parent itself
  differs from the saved cluster checkpoint by at most 3.553e-15.
- Population reconstruction L1 = 4.297e-15; birth gap =
  2.359e-16; no feasibility-projection mass.
- Zero occupied budget violations; maximum occupied budget excess 1.776e-15.
- Zero occupied downward value steps; probability arrays finite and in [0,1].
- Complete bounded run: 139.3 seconds, including checking
  and saving; full two-shock solve 111.0 seconds.

The existing local runtime could not load the NumPy2 checkpoint, so an isolated
Python3.12/NumPy2.2.6/SciPy1.15.3/Numba0.61.2 environment was used. Cluster login
was rejected; no cluster job was submitted. Initial failed attempts remain
preserved. Final receipts are only under `final_comparison/`.

No figures were generated, following the author's explicit preference.
Reproduction instructions: `docs/model/e5f_two_shock_test.md` on the experimental
branch. Remaining work: evaluate the unchanged twelve empirical moments, clear
markets and recalibrate this specification before making policy claims.
