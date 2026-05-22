# Income Risk / Mortgage Account Smoke Report

Verdict: **yellow**

## What Ran

The copied benchmark solver ran once at the current direct benchmark theta and
then constructed finite diagnostic grids for \(z\) and \(\mu\). The grids are
non-degenerate, but the Bellman recursion and forward distribution in this
overnight prototype do **not** yet optimize over those two states. This is a
diagnostic state augmentation, not an accepted structural HANK/mortgage branch.

## Solve Status

- accepted equilibrium from copied baseline solver: `True`
- convergence reason: `strict_tol`
- best equilibrium error: `0.00048043163102334807`
- elapsed seconds: `47.49`
- runtime category: `moderate` for one copied benchmark solve; full structural
  \(z,\mu\) state expansion remains `project-scale`
- SMM loss on this smoke solve: `23.848`

## State Sizes

| Dimension | Count |
|---|---:|
| liquid wealth b | 80 |
| tenure d | 7 |
| locations | 2 |
| ages | 60 |
| completed-fertility states | 4 |
| child-age states | 7 |
| earnings z states | 3 |
| mortgage-account mu states | 5 |

## Branch Diagnostics

- non-degenerate z states used in diagnostic distribution: `3`
- non-degenerate mu states used in diagnostic distribution: `5`
- maximum diagnostic default hazard: `0.000415457`
- average bad-credit share at age 35: `0.0723944`
- origination shares by credit status: `{"G": 0.9199999999999997, "B": 0.07999999999999993}`

## Interpretation

This fails the full Branch 1 acceptance criterion because the new states are
not yet inside the Bellman equation. It is still useful as a smoke test of the
state definitions, diagnostic accounting, and report pipeline. The next
implementation step is to expand value, policy, and distribution arrays from
`(b,d,i,a,n,s)` to `(b,d,i,a,n,s,z,mu)` and to feed \(y_a(z)\), mortgage
payments, amortization, default, and credit recovery directly into the copied
solver.
