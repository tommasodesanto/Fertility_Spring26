# Income Risk / Mortgage Account Prototype

This folder is an isolated overnight branch test. It contains a copied
`dt_cp_model/` package and does not modify the live implementation under
`code/model/dt_cp_model`.

The smoke driver is:

```bash
cd overnight_variants/2026-05-22_income_mortgage_risk
/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/.venv/bin/python run_income_mortgage_risk_smoke.py --quiet
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
