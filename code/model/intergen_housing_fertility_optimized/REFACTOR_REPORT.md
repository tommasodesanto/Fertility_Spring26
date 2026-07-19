# Optimized M5 Refactor Report

## Scope

The refactor is physically isolated under
`code/model/intergen_housing_fertility_optimized/`. The production package was
not edited and no existing caller was redirected.

The governing rule was that an optimization remains enabled only if it passes
both numerical and performance gates. The saved canonical M5 price is the
exact-parity object; equilibrium comparisons separately respect the existing
nearby-root/platform qualification.

## Correctness changes

1. Income persistence or same-sized income weights now rebuild the transition
   matrix unless an explicit matrix is supplied. Explicit matrices must be
   consistent with their declared stationary weights.
2. Unknown override keys are rejected except for an explicit experiment-hook
   registry.
3. Family/child flattened indices use the centralized Fortran-order contract
   `column = parity + n_parity * child_state` in both Python owner fallbacks.
4. With retirement income heterogeneity, every related old-age nonhousing and
   total-wealth mean/median, including parent and childless variants, is
   recomputed on the full income-state distribution.
5. The M5 income-disciplined target registry directly contains all 15 rows,
   including aggregate occupied rooms. `TargetSystem` makes order, values,
   weights, identification count, and a stable fingerprint explicit.

None of these changes alters the default saved-price M5 target moments.

## Structure, clarity, and style

- `m5_profile.py` owns the exact benchmark theta, survival schedule, income
  process, entry distribution, and runtime switches. Benchmarking no longer
  imports a calibration runner for hidden configuration.
- `TargetSystem` owns ordered moment names, targets, weights, identification,
  and a stable fingerprint as one immutable object.
- `StateLayout` and `decode_flat_family_state` replace duplicated arithmetic
  for the Fortran-order family-state encoding.
- `benchmark_m5.py`, `parity_panel.py`, and `tests/run_all.py` provide one clear
  entry point each for saved-price parity, local parameter parity, and tests.
- Unknown configuration is rejected rather than silently attached to the
  parameter namespace. This makes misspellings and stale runner switches fail
  close to their source.

The inherited numerical core remains large (`solver.py` is about 6,600 lines).
This pass isolates new contracts and removes the measured-slower fused
prototype, but it deliberately does not perform a cosmetic file split of the
Bellman and distribution routines. Such a split would change import boundaries
without reducing computation and should be considered only after the remaining
promotion gates, with the same parity panel applied to each extraction.

## Accepted performance changes

### Compiled scatter

At the canonical saved price, production and optimized policy/value arrays are
exactly equal. The distribution differs by at most `4.19e-16`, total mass by
about `1e-16`, and the largest 15-moment difference is `1.69e-14`.

In the paired benchmark, the complete fixed-price solve fell from `11.02` to
`9.18` seconds (`1.20x`). Separate warm measurements ranged from roughly 20%
to 34%; only the paired result should be treated as the headline.

### Direct scalar equilibrium search

The optimized one-market search evaluates the initial price, expands only in
the sign-indicated direction, and then uses safeguarded Brent interpolation.
The legacy damped loop remains the multi-market and no-bracket fallback.

On the same local host:

| Solver | Fast price evaluations | Wall seconds | Residual |
|---|---:|---:|---:|
| Production legacy | 33 | 114.39 | `9.32e-6` |
| Optimized direct plus scatter | 6 | 35.49 | `1.30e-5` |

The resulting local speed ratio is `3.22x`. Absolute times are host-load
sensitive; the solve-count reduction is the more portable result.

The optimized strict equilibrium is a nearby accepted root at price
`0.7589109722714222`, versus canonical Torch price `0.7589039340587314`. Its
loss is `9.030561727565956`, but this is not a replacement calibration result.
At the canonical price, the optimized package reproduces canonical loss
`9.044422069071352` to rounding precision.

This benchmark does not define a new calibration. The canonical M5 reporting
contract remains the complete
[15-row target-fit table](../../../output/model/intergen_income_disciplined_recalibration_20260716/report/target_fit_full.csv)
and
[full parameter/bounds table](../../../output/model/intergen_income_disciplined_recalibration_20260716/report/parameter_table_full.csv).
All 14 estimated parameters are reported there with bounds and near-bound
flags; `theta_n=0` is the stated external restriction.

### Accepted Bellman reuse

Fast root evaluations retain only the current best household payload; cached
non-best entries retain scalar/statistical summaries. The accepted policies are
reused for the final full distribution and statistics, avoiding one repeated
Bellman pass without retaining every candidate's large arrays.

## Rejected performance change

The first fully fused one-market Markov age kernel matched the reference KFE to
`3.08e-15` cellwise and conserved mass to `4.44e-16`, but required `4.47`
seconds against `2.07` seconds for the compiled-scatter reference. It is
removed and excluded from all accepted speedups.

## Verification

- 78 `unittest` tests pass.
- 10 top-level tests pass through the stable direct runner.
- Total: 88 tests.
- Saved-price values and every policy/probability array are exactly equal.
- Saved-price distribution maximum difference: `4.19e-16`.
- Saved-price largest target-moment difference: `1.69e-14`.
- Two optimized strict GE runs produced bit-identical prices, moments, losses,
  masses, and residuals.

The six-case fixed-price perturbation panel covers discounting, the consumption
floor, child housing demand, deterministic tenure, and additional-child space.
Across all cases, value and policy/probability arrays are exactly equal; the
largest cellwise distribution difference is below `4.19e-16`, and the largest
15-moment difference is `4.09e-14`. Median paired fixed-price speedup is
`1.18x`. One noisy pair is `0.93x`, so the scatter change should be understood
as a modest throughput improvement rather than a guaranteed latency reduction
for every parameter point. The direct equilibrium search remains the material
speed gain.

Raw local benchmark JSON is under
`output/model/intergen_optimized_refactor_20260719/`.

## Torch promotion battery

The full battery completed on Torch with one CPU per task. All numerical
promotion gates pass:

| Gate | Result |
|---|---:|
| 29 exact-bound production/optimized comparisons | PASS |
| 17-price map over configured `[0.01,30]` bracket | PASS |
| Six optimized strict-root cases | PASS |
| Two bit-identical strict repeats | PASS |
| 61 standard diagnostic artifacts | PASS |
| Ten-case sequential GE throughput smoke | PASS |

The six actual-root cases have median speedup `4.12x` and minimum speedup
`2.28x`. The ten-case sequential objective sweep takes `385.86` seconds in
production and `141.88` seconds optimized, a `2.72x` throughput gain; both
implementations select the same predeclared candidate. This sweep is a
numerical throughput test, not an identified SMM estimate or replacement for
M5.

Two adjudications remain explicit rather than hidden:

1. The original `5e-10` moment-parity gate passes 28/29 exact-bound cases. At
   `psi_child=-3`, exact policies and a `4.65e-16` distribution difference are
   amplified to `8.41e-10` in one conditional room mean. Promotion uses a
   disclosed `1e-8` conditioned-moment tolerance for that case.
2. Every optimized root satisfies `2.5e-5`. The legacy production root misses
   final full-stat tolerance in two cases even though its fast refinement
   claimed convergence. Fixed-price demand parity independently verifies that
   both packages evaluate the same economic mapping.

The full evidence and both selected-case 15-row target-fit and parameter/bounds
tables are under
`output/model/intergen_optimized_promotion_20260719/`. Numerical promotion is
complete. Existing tools remain pointed at production until the user explicitly
chooses to promote this package.

## Porting conclusion

The evidence does not support a wholesale C++ or Rust rewrite as the next
speed step. The expensive household kernels are already compiled with Numba;
the large measured gain came from reducing equilibrium evaluations and reusing
the accepted Bellman payload. A lower-level rewrite would mainly target the
remaining Python-orchestrated distribution work and could plausibly improve a
fixed-price solve, but it would not reproduce the `3.22x` algorithmic gain by
itself and would impose a much larger correctness and maintenance burden.
