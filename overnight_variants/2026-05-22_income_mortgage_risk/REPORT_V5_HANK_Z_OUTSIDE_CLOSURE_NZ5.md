# Income Risk V5 HANK-z Outside-Option Closure Report

Verdict: **yellow**

## What Changed

This run keeps the structural HANK earnings state \(z\) and switches the copied
GE loop from the renewal-valve closure to the paper-facing outside-option scale
closure. The closure implements
\[
S E_0(p)=q^E(p)\left[M+S B_0(p)\right],
\qquad
S(p)=\frac{q^E(p)M}{E_0(p)-q^E(p)B_0(p)}.
\]
For the benchmark steady state, \(S=1\) is imposed inside the GE loop. At each
candidate composition, the script computes \(E_0(p)\) and \(B_0(p)\),
calibrates \(\bar W^E\) to the target \(q^E\) at current entry values, and sets
\[
M = E_0(p)/q^E(p)-B_0(p).
\]
This replaces the earlier outer normalization pass. The reported final \(M\)
and \(\bar W^E\) are the objects to hold fixed in counterfactuals.

## Benchmark Normalization

- target \(q^E\): `0.9000`
- calibrated outside value \(\bar W^E\): `-1.50414e+06`
- entry logit scale \(\kappa_E\): `1e+06`
- residual outside-born flow \(M\): `0.00580643`
- entry flow per unit scale \(E_0\): `0.0166667`
- mature city-born flow per unit scale \(B_0\): `0.0127121`
- accounting scale factor \(S\): `1`
- scale residual: `0.000e+00`

## Final GE Status

- closure: `outside_option_benchmark_normalized`
- accepted: `True`
- strict converged: `True`
- convergence reason: `strict_tol`
- iterations completed: `12`
- best equilibrium error: `0.000346432`
- final equilibrium error: `0.000346432`
- prices: `[0.511963296326644, 0.6124693819761028]`
- \(z\) states: `5`
- \(b\) states: `30`
- final scale factor \(S\): `1`
- final city-entry probability \(q^E\): `0.9`
- final outside probability: `0.1`
- finite scale: `True`
- elapsed seconds: `83.10`
- SMM loss against live targets: `311.3`

## Normalization Passes

No outer normalization pass was used. The benchmark normalization is imposed
directly inside each GE iteration.

## Moment Table

| Moment | Target | Benchmark | V5 outside closure |
|---|---:|---:|---:|
| `tfr` | 1.700 | 1.898 | 1.558 |
| `childless_rate` | 0.150 | 0.145 | 0.351 |
| `mean_age_first_birth` | 26.000 | 33.535 | 34.793 |
| `tfr_gradient` | 0.133 | 0.119 | -0.239 |
| `own_rate` | 0.627 | 0.643 | 0.575 |
| `own_gradient` | 0.170 | 0.139 | -0.082 |
| `own_family_gap` | 0.110 | 0.114 | 0.383 |
| `prime_childless_renter_median_rooms` | 4.000 | 6.365 | 5.767 |
| `prime_childless_owner_median_rooms` | 6.000 | 6.800 | 6.800 |
| `housing_increment_0to1` | 0.664 | 0.441 | 0.591 |
| `housing_increment_1to2` | 0.566 | 0.192 | 1.492 |
| `young_liquid_wealth_to_income` | 0.600 | 0.527 | 1.215 |
| `center_share_nonparents` | 0.494 | 0.405 | 0.316 |
| `center_share_newparents` | 0.416 | 0.382 | 0.393 |
| `migration_rate` | 0.032 | 0.035 | 0.034 |
| `old_age_own_rate` | 0.863 | 0.947 | 0.709 |
| `old_age_parent_childless_gap` | 0.070 | 0.062 | 0.291 |
| `inv_pop_share_C` | 0.450 | 0.441 | 0.409 |
| `inv_rent_ratio_C_over_P` | 1.140 | 1.182 | 1.196 |

## Read

Auditing against `latex/model_writeup.tex` and the copied solver logic, the
outer re-solve was unnecessarily complicated for the benchmark. The writeup
normalizes the baseline at \(S=1\) and sets \(M\) residually from the baseline
objects. Numerically, the cleaner benchmark implementation is to impose that
residual normalization at the current candidate equilibrium/composition and use
normalized housing demand. A fixed-\(M\) scale formula is the counterfactual
object, not the object that should move the benchmark away from its own
normalization.

The run is still yellow because the un-recalibrated economics are not
acceptable yet: fertility is too low and too late, geography and ownership
gradients flip sign, and liquid wealth is too high. The mortgage account state
is also still absent from this HANK-z-only branch.

## Closure Lessons

The entry margin needs its own scale \(\kappa_E\). Using the incumbent
within-city location scale \(\kappa_\ell\) made the outside probability jump to
zero or one because entry values are lifetime-utility objects with very large
levels. The accepted run therefore uses an explicitly supplied
\(\kappa_E=1e+06\). The Rouwenhorst grid also requires the
true stationary Markov distribution; this script now reports the binomial
stationary weights for \(N_z=5\), not uniform smoke weights.
