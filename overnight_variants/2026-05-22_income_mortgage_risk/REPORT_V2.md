# Income Risk / Mortgage Account V2 Scenario Report

Verdict: **yellow**

## What Changed Relative To V1

V1 only overlaid \(z\) and \(\mu\) after a benchmark solve. V2 solves a finite
set of partial-equilibrium scenarios where the new objects enter choices:

- \(z\) scales working-age earnings through \(y_a(z)\);
- good and bad credit change the financed share \(\phi_g\);
- active mortgage accounts add a compact owner-payment wedge equal to
  \(\rho_g m\) as an owner user-cost increment.

This is still not a full branch. There is no joint Bellman state
\((b,d,i,a,n,s,z,\mu)\), no endogenous default choice inside the period, and no
single distribution transition over \(z\) and \(\mu\). It is a scenario-mixture
approximation designed to test whether those margins move choices enough to
justify the real state expansion.

## Scenario Grid

- \(z\) states: `[-0.28, 0.0, 0.28]`
- \(\mu\) states: `['active_high_good', 'active_low_bad', 'active_low_good', 'no_mortgage_bad', 'no_mortgage_good']`
- scenario solves: `15`
- elapsed seconds: `44.76`
- runtime category: `moderate` for the scenario mixture; full structural
  expansion remains `project-scale`
- scenario-mixture SMM loss: `80.0675`

## Movement

| Object | Value |
|---|---:|
| ownership spread across \(z\) | 0.073 |
| TFR spread across \(z\) | 2.487 |
| max default-rate proxy | 0.973 |
| aggregate TFR | 1.532 |
| aggregate ownership | 0.531 |
| aggregate young liquid wealth / income | 0.443 |

## Read

The scenario mixture produces meaningful variation in ownership and fertility
across \(z\) and \(\mu\), which is an improvement over V1. It is not acceptable
as a branch result: aggregate childlessness, first-birth timing, geography, and
old-age parent-childless ownership deteriorate. The yellow verdict means only
that the risk/account margins have enough bite to justify a real structural
prototype. That next implementation should expand value, policy, and
distribution arrays over \(z\) and \(\mu\), then add the default decision before
fertility/location/tenure choices.
