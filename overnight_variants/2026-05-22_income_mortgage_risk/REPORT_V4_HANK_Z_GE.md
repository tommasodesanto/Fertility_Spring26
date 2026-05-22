# Income Risk V4 HANK-z Full Equilibrium Report

Verdict: **yellow**

## What This Is

This is the full-equilibrium correction to V3. The copied branch solves prices
and entry shares with the structural HANK earnings state \(z\) inside the
household problem and forward distribution. It is no longer a fixed-price PE
test. The household state is \((b,d,i,a,n,s,z)\), with
\(y_{ia}(z)=(1-\tau_{pay})w_i e_a\exp(z)\) for working ages and a common
retirement mapping.

The mortgage-account state \(\mu\) is still not structural. This run answers
the narrower question: can the classic HANK \(z\)-state version clear the
copied model's equilibrium loop on a coarse grid?

## Equilibrium Status

- accepted: `True`
- strict converged: `True`
- convergence reason: `strict_tol`
- iterations completed: `15`
- best equilibrium error: `0.000167292`
- best iteration: `15`
- final equilibrium error: `0.000167292`
- prices: `[0.505855651256028, 0.6001949358069548]`
- \(z\) states: `3`
- \(b\) states: `30`
- elapsed seconds: `85.56`
- runtime category: `expensive`
- SMM loss against live targets: `45.9235`

## HANK Diagnostics

| Object | Value |
|---|---:|
| ownership spread across \(z\) | 0.109 |
| fertility-choice spread across \(z\) | 0.034 |
| aggregate TFR | 1.858 |
| aggregate ownership | 0.639 |
| young liquid wealth / income | 0.935 |

## Moment Table

| Moment | Target | Benchmark | V4 GE |
|---|---:|---:|---:|
| `tfr` | 1.700 | 1.898 | 1.858 |
| `childless_rate` | 0.150 | 0.145 | 0.181 |
| `mean_age_first_birth` | 26.000 | 33.535 | 34.318 |
| `tfr_gradient` | 0.133 | 0.119 | 0.087 |
| `own_rate` | 0.627 | 0.643 | 0.639 |
| `own_gradient` | 0.170 | 0.139 | -0.066 |
| `own_family_gap` | 0.110 | 0.114 | 0.154 |
| `prime_childless_renter_median_rooms` | 4.000 | 6.365 | 6.180 |
| `prime_childless_owner_median_rooms` | 6.000 | 6.800 | 6.800 |
| `housing_increment_0to1` | 0.664 | 0.441 | 0.415 |
| `housing_increment_1to2` | 0.566 | 0.192 | 0.402 |
| `young_liquid_wealth_to_income` | 0.600 | 0.527 | 0.935 |
| `center_share_nonparents` | 0.494 | 0.405 | 0.375 |
| `center_share_newparents` | 0.416 | 0.382 | 0.372 |
| `migration_rate` | 0.032 | 0.035 | 0.035 |
| `old_age_own_rate` | 0.863 | 0.947 | 0.772 |
| `old_age_parent_childless_gap` | 0.070 | 0.062 | 0.004 |
| `inv_pop_share_C` | 0.450 | 0.441 | 0.433 |
| `inv_rent_ratio_C_over_P` | 1.140 | 1.182 | 1.186 |

## Read

This is the equilibrium object the branch needed: \(z\) is a real Markov state
and the price loop moves with the \(z\)-state distribution. The result is still
yellow, not green, because the coarse grid is not recalibrated and the mortgage
account is absent. The next live implementation should add only this HANK-z
core first, then decide separately whether a compact account state is worth the
extra state-space cost.
