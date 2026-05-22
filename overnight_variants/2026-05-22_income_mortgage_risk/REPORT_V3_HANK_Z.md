# Income Risk V3 HANK-z Structural Report

Verdict: **yellow**

## What Changed

This pass implements the classic HANK object the earlier passes were missing:
a finite idiosyncratic earnings state \(z\) with a Markov transition matrix
\(\Pi_z\). The copied solver now carries value functions, policies, fertility
probabilities, location probabilities, tenure choices, and the forward
distribution over \((b,d,i,a,n,s,z)\). Continuation values average over
\(\Pi_z\) after child aging. Working-age income is
\(y_{ia}(z)=(1-\tau_{pay})w_i e_a\exp(z)\); retirement uses the copied
common pension mapping.

The mortgage account \(\mu\) is not structural in V3. That is intentional for
this correction: the HANK earnings state is now real, while defaultable account
dynamics remain the next layer because \(\mu'\) must depend on purchase, sale,
amortization, and default timing.

## Status

- solve: fixed-price partial equilibrium at copied benchmark prices
- \(z\) states: `3`
- \(z\) grid: `[-0.28, 0.0, 0.28]`
- stationary \(z\) distribution: `[0.3333333333333333, 0.3333333333333333, 0.3333333333333333]`
- \(b\) states: `30`
- tenure states: `7`
- locations: `2`
- ages: `60`
- fertility states: `4`
- child-age states: `7`
- elapsed seconds: `12.62`
- runtime category: `moderate`
- SMM loss against live targets: `48.0415`

## HANK Diagnostics

| Object | Value |
|---|---:|
| ownership spread across \(z\) | 0.109 |
| fertility-choice spread across \(z\) | 0.035 |
| aggregate TFR | 1.880 |
| aggregate ownership | 0.645 |
| young liquid wealth / income | 0.929 |

## Moment Table

| Moment | Target | Benchmark | V3 |
|---|---:|---:|---:|
| `tfr` | 1.700 | 1.898 | 1.880 |
| `childless_rate` | 0.150 | 0.145 | 0.175 |
| `mean_age_first_birth` | 26.000 | 33.535 | 34.257 |
| `tfr_gradient` | 0.133 | 0.119 | 0.085 |
| `own_rate` | 0.627 | 0.643 | 0.645 |
| `own_gradient` | 0.170 | 0.139 | -0.079 |
| `own_family_gap` | 0.110 | 0.114 | 0.160 |
| `prime_childless_renter_median_rooms` | 4.000 | 6.365 | 6.238 |
| `prime_childless_owner_median_rooms` | 6.000 | 6.800 | 6.800 |
| `housing_increment_0to1` | 0.664 | 0.441 | 0.412 |
| `housing_increment_1to2` | 0.566 | 0.192 | 0.401 |
| `young_liquid_wealth_to_income` | 0.600 | 0.527 | 0.929 |
| `center_share_nonparents` | 0.494 | 0.405 | 0.393 |
| `center_share_newparents` | 0.416 | 0.382 | 0.388 |
| `migration_rate` | 0.032 | 0.035 | 0.035 |
| `old_age_own_rate` | 0.863 | 0.947 | 0.803 |
| `old_age_parent_childless_gap` | 0.070 | 0.062 | 0.006 |
| `inv_pop_share_C` | 0.450 | 0.441 | 0.446 |
| `inv_rent_ratio_C_over_P` | 1.140 | 1.182 | 1.182 |

## Read

This is the correct computational direction for Branch 1: idiosyncratic labor
income risk is now an actual discrete state, not a post-solve label or a set of
separate scenario economies. The test should still be read as yellow rather
than green because it is fixed-price PE and does not yet include the structural
mortgage-account state \(\mu\). The next Branch 1 pass should add \(\mu\) as a
second finite state with tenure-dependent account transitions, then rerun this
same target table.
