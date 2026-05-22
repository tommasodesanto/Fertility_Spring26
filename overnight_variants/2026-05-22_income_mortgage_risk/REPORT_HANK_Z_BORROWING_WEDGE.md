# HANK-z Borrowing-Wedge Diagnostic

Verdict: **yellow**

## What This Tests

This keeps Branch 1 as one additional state, \(z\), and does not add a
separate mortgage/default account. The question is whether the existing
housing-finance wedge,
\[
(1-\phi)p_iH_k \quad \text{and} \quad b' \ge -\phi p_iH_k,
\]
is used once households face persistent earnings risk.

## Solve Status

- GE accepted: `True`
- convergence reason: `strict_tol`
- final equilibrium error: `0.000167292`
- prices: `[0.505855651256028, 0.6001949358069548]`
- \(z\) grid: `[-0.28, 0.0, 0.28]`
- \(b\) states: `30`
- elapsed seconds: `108.56`
- runtime category: `expensive`
- SMM loss against live targets: `45.9235`

## Wedge Diagnostics

| Diagnostic | Low \(z\) | High \(z\) |
|---|---:|---:|
| ownership rate | 0.521 | 0.630 |
| renter-to-owner purchase share | 0.015 | 0.091 |
| renter feasible for starter down payment | 0.336 | 0.590 |
| renter mass near starter down payment | 0.881 | 0.724 |
| purchase mean down-payment slack | 0.982 | 0.920 |
| owner \(b'\) near borrowing floor | 0.089 | 0.029 |

## Moment Table

| Moment | Target | Benchmark | HANK-z GE |
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

The existing \(\phi\)-based wedge is economically active. High-\(z\) renters
are much more likely to buy than low-\(z\) renters, and they are much more
likely to be feasible for the starter down payment. Low-\(z\) owners also have
more mass close to the borrowing floor. That is the validation we wanted: the
current model already has a mortgage-like collateral/liquidity wedge.

The diagnostic is not green for the full branch. Actual purchasers are not
tightly bunched at the down-payment boundary, young liquid wealth remains too
high, and the ownership gradient is still wrong-signed. So this validates the
direction of using HANK-\(z\) with the existing \(\phi\) wedge, but it points to
recalibration and grid/transition discipline rather than adding a rich
mortgage/default state.
