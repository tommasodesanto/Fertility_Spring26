# Moment-matched period proxy

A mean-one lognormal process can closely approximate the level covariances of four-year average earnings while retaining a persistent component and an independent transitory component. This is an approximation check, not an adopted income process or a household-model result.

Writing $C_k$ for the level covariance at lag $k$, match the exact block-average covariances at lags zero, one and two:

\[
\rho=\frac{\log(1+C_2)}{\log(1+C_1)},\qquad
V_p=\frac{\log(1+C_1)}{\rho},\qquad
V_e=\log(1+C_0)-V_p.
\]

The implied period persistence is $\rho=0.8863302844$, persistent log variance $V_p=0.6941543660$ and transitory log variance $V_e=0.0737459191$. All are admissible. The implied persistent innovation variance is $V_p(1-\rho^2)$; $V_p$ itself is the stationary persistent variance.

| Lag (four-year periods) | Exact block covariance | Matched proxy | Existing continuous endpoint proxy |
|---|---:|---:|---:|
| 0 | 1.15523612 | 1.15523612 | 1.20489480 |
| 1 | 0.85011914 | 0.85011914 | 0.84803365 |
| 2 | 0.72515127 | 0.72515127 | 0.72350210 |
| 4 | 0.53480029 | 0.53478315 | 0.53373465 |

The lag-four error is $-0.0032\%$. This comparison precedes finite-grid discretization: using 15 discrete states introduces a separate approximation. Matching three covariances does not identify shock timing, tails or conditional distributions, and therefore does not establish equivalence for household choices. The comparison also conditions on the existing annual-process estimates, whose empirical measurement contract remains a separate decision.

The lead checked the receipt hash and recomputed every proxy covariance from the displayed formulas. Full precision, input hashes and relative errors are in [moment_matched_period_proxy.json](moment_matched_period_proxy.json). No household solve, target or calibration change occurred.
