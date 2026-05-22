# Income Risk / Mortgage Account Prototype

This folder is an isolated overnight branch test. It contains a copied
`dt_cp_model/` package and does not modify the live implementation under
`code/model/dt_cp_model`.

The smoke driver is:

```bash
cd overnight_variants/2026-05-22_income_mortgage_risk
/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/.venv/bin/python run_income_mortgage_risk_smoke.py --quiet
```

The second-pass scenario driver is:

```bash
cd overnight_variants/2026-05-22_income_mortgage_risk
/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/.venv/bin/python run_income_mortgage_risk_v2_scenarios.py --quiet --nb 40
```

The third-pass HANK earnings-risk driver is:

```bash
cd overnight_variants/2026-05-22_income_mortgage_risk
/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/.venv/bin/python run_income_mortgage_risk_v3_hank_z.py --quiet --nb 30 --nz 3
```

The full-equilibrium HANK earnings-risk driver is:

```bash
cd overnight_variants/2026-05-22_income_mortgage_risk
/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/.venv/bin/python run_income_mortgage_risk_v4_hank_z_ge.py --quiet --nb 30 --nz 3 --max-iter-eq 35
```

## Implemented Objects

- Earnings grid \(z\in\{-0.28,0,0.28\}\) with a symmetric persistent Markov
  matrix.
- Mortgage-account grid \(\mu\) with five encoded states:
  no mortgage/good credit, no mortgage/bad credit, active low/good, active
  high/good, and active low/bad.
- Diagnostic mortgage objects: amortization rate, good/bad payment rates,
  good/bad financed shares, default penalty, and bad-credit recovery.
- A copied-solver helper that attaches an augmented diagnostic distribution
  over \((z,\mu)\) after one benchmark equilibrium solve.

## Important Limitation

The overnight implementation does not yet put \(z\) and \(\mu\) inside the
Bellman recursion. The household policies are still solved on the baseline
state \(x=(b,d,i,a,n,s)\). The augmented objects are diagnostic smoke-test
accounting and should not be treated as an accepted structural HANK/mortgage
implementation.

The generated `REPORT.md` states this explicitly and gives the branch verdict.

## Second-Pass V2

`run_income_mortgage_risk_v2_scenarios.py` improves on V1 by making the new
objects affect solved choices in a finite set of partial-equilibrium scenarios.
The script solves 15 economies: 3 earnings states times 5 mortgage-account
states. Earnings scale \(y_a(z)\), bad credit lowers the financed share
\(\phi_g\), and active mortgage accounts add a compact owner-payment wedge
\(\rho_g m\).

This still is not a full structural implementation because there is no single
Bellman recursion over \((b,d,i,a,n,s,z,\mu)\). See `REPORT_V2.md`,
`results_income_mortgage_risk_v2.csv`, and
`diagnostics_income_mortgage_risk_v2.csv`.

## Third-Pass V3 HANK-z

`run_income_mortgage_risk_v3_hank_z.py` implements the standard HANK-style
idiosyncratic earnings process as a true state. The copied solver carries
value functions, policies, fertility probabilities, location probabilities,
tenure choices, and the forward distribution over \((b,d,i,a,n,s,z)\). The
earnings state follows the finite Markov chain reported in
`diagnostics_income_mortgage_risk_v3_hank_z.csv`, and continuation values
average over \(\Pi_z\) after child aging.

This is the correction to the earlier Branch 1 mistake: \(z\) is no longer a
post-solve label or a scenario loop. The output remains a partial-equilibrium
prototype at copied benchmark prices, and \(\mu\) is not yet structural. Adding
\(\mu\) requires tenure-dependent account transitions for purchase,
amortization, sale, and default. See `REPORT_V3_HANK_Z.md`,
`results_income_mortgage_risk_v3_hank_z.csv`, and
`diagnostics_income_mortgage_risk_v3_hank_z.csv`.

## Fourth-Pass V4 HANK-z Full Equilibrium

`run_income_mortgage_risk_v4_hank_z_ge.py` runs the same structural \(z\)-state
household problem inside the copied equilibrium price/entry loop. This
supersedes the V3 fixed-price result as the relevant Branch 1 smoke evidence.
On the coarse `Nb=30`, `Nz=3` run, it accepted at strict tolerance in 15
iterations and took 85.56 seconds. See `REPORT_V4_HANK_Z_GE.md`,
`results_income_mortgage_risk_v4_hank_z_ge.csv`,
`diagnostics_income_mortgage_risk_v4_hank_z_ge.csv`, and
`diagnostics_income_mortgage_risk_v4_hank_z_ge_trace.csv`.
