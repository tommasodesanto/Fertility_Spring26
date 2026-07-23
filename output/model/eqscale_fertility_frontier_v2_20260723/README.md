# Fertility-frontier v2: imposed equivalence scales under FL x HSV income (2026-07-23)

Torch array `14679019` (smoke `14678986`), 147/147 cells completed, max GE
residual `1.2e-04`. Design: {power (SSK), sqrt, linear_g05} x 7 psi_child x
7 kappa_fert at the certified E3 winner, L4 literal-parity conventions, and
the Floden-Linde x HSV after-tax income process (rho 0.9136, annual
innovation s.d. 0.1690). The linear_g05 layer isolates the income-process
change; the v1 baseline (linear scale, old income pair) is
`output/model/eqscale_fertility_frontier_20260722/`.

Pre-registered question: does the concave imposed scale unlock completed
fertility 1.918 AND childlessness 0.188 jointly?

**Answer: no.** Zero cells within 10 percent of both stocks under any form.
Best joint-distance cells:

| Layer | Best cell | psi | kappa_f | CF | childless | joint dist |
|---|---|---:|---:|---:|---:|---:|
| v1 baseline (linear, old income) | 107 | 0.0 | 1.0 | 1.262 | 0.173 | 0.1233 |
| power (SSK) | 4 | -0.5 | 3.0 | 1.254 | 0.173 | 0.1261 |
| sqrt | 53 | -0.5 | 3.0 | 1.258 | 0.172 | 0.1255 |
| linear_g05 (income control) | 107 | 0.0 | 1.0 | 1.331 | 0.148 | 0.1388 |

The frontier is invariant to scale curvature (linear/power/sqrt), scale
generosity (gamma_e 0.0085 to 0.5), and the income-process change. At every
best joint cell the 3+ share is 0.07-0.10 while hitting CF 1.918 at
childlessness 0.188 requires parents to average 2.36 children. The missing
margin is continuation intensity among parents holding entry fixed; a
single psi with iid per-period taste noise ties the two margins in every
configuration tested. This is the pre-registered evidence that the binding
failure is the preference structure itself, not the needs-scale form or the
income process.
