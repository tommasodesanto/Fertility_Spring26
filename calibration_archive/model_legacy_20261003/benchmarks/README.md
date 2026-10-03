# Benchmark Notes

Generated on 2026-04-25.

## Reduced Fast GE, One Full Iteration

Python command:

```bash
.venv/bin/python -m dt_cp_model.cli benchmark --setup fast --max-iter-eq 1 --force-full --quiet --json benchmarks/fast_ge_iter1.json
```

MATLAB command:

```bash
/Applications/MATLAB_R2025b.app/bin/matlab -batch "addpath('code/model/tools'); matlab_reference_fast_iter1"
```

Comparison:

```bash
.venv/bin/python tools/compare_benchmarks.py benchmarks/fast_ge_iter1.json benchmarks/matlab_fast_ge_iter1.json
```

Key results:

- Python total wall time: `15.81s`
- MATLAB total wall time: `43.12s`
- Python full Bellman time: `6.22s`
- MATLAB full Bellman time: `29.34s`
- Max absolute moment difference in the comparison table: `0.00263`
- Prices match exactly in this reduced benchmark.
- Distribution mass is conserved to floating-point precision.

This benchmark is not a converged equilibrium. It is a parity and speed check
for the translated Bellman, transition, and moment machinery at live fast-setup
dimensions.

## Reduced Fast GE, Live Stopping Rule

Python command:

```bash
.venv/bin/python -m dt_cp_model.cli benchmark --setup fast --max-iter-eq 120 --json benchmarks/python_fast_ge_converged.json
```

MATLAB command:

```bash
/Applications/MATLAB_R2025b.app/bin/matlab -batch "addpath('code/model/tools'); matlab_reference_fast_converged"
```

Both solvers followed the same convergence path and hit the live soft-accept
criterion at iteration `26`:

- Python `best_eq_error`: `1.3249e-3`
- MATLAB printed `best_err=1.33e-3`
- live `fast` tolerance: `1e-3`

This is the MATLAB solver's existing "accepted equilibrium" rule, not a strict
`err < tol` termination.

Key results:

- Python total wall time: `151.23s`
- MATLAB total wall time: `182.95s`
- Python full Bellman time: `29.48s` over `9` calls
- MATLAB full Bellman time: `126.85s` over `9` calls
- Python distribution time: `85.60s` over `26` calls
- MATLAB distribution time: `29.87s` over `26` calls
- Max absolute moment difference in the comparison table: `0.00531`
- Prices differ by less than `5e-6`.
- Distribution mass is conserved to floating-point precision.

Interpretation:

- The Python port matches MATLAB closely on a converged-to-live-rule fast GE
  benchmark.
- The Bellman block is already faster in Python with Numba kernels.
- The distribution block is now the main Python bottleneck and should be the
  first target for compiled scatter-add kernels.

## Direct Theta X0, Fast Setup

This is the normal model boundary for calibration work:

```bash
.venv/bin/python -m dt_cp_model.cli solve-theta --setup fast --theta x0 --max-iter-eq 120 --quiet --json benchmarks/solve_theta_x0_fast_compiled_forward.json
```

MATLAB reference:

```bash
/Applications/MATLAB_R2025b.app/bin/matlab -batch "addpath('code/model/tools'); matlab_reference_solve_theta_x0_fast"
```

Comparison:

```bash
.venv/bin/python tools/compare_benchmarks.py benchmarks/solve_theta_x0_fast_compiled_forward.json benchmarks/matlab_solve_theta_x0_fast.json
```

Key results:

- Python total wall time: `196.57s`
- MATLAB total wall time: `268.83s`
- Python full Bellman time: `75.36s` over `8` calls
- MATLAB full Bellman time: `187.73s` over `8` calls
- Python Howard eval time: `82.42s` over `18` calls
- MATLAB Howard eval time: `32.09s` over `18` calls
- Python distribution time: `30.14s` over `26` calls
- MATLAB distribution time: `41.52s` over `26` calls
- Max absolute moment difference in the comparison table: `0.00547`
- Prices differ by less than `7e-5`.
- Distribution mass is conserved to floating-point precision.

Interpretation:

- The direct `theta -> solve_theta -> moments` path works and is now faster
  than the matching MATLAB benchmark for the fast accepted-equilibrium case.
- Python wins in the full Bellman solve and now also in the repeated
  fast-stat distribution propagation.
- MATLAB is still faster in Howard eval, so cached/policy-evaluation work is
  the next speed target before larger SMM runs.

## Distribution Optimization Pass 1

Changes:

- Intermediate equilibrium iterations now skip event-study housing diagnostics
  when `fast_stats=true`; these diagnostics do not feed price or entry updates
  and are still computed in the final full-stat pass.
- Default redistribution changed from `np.add.at` to `np.bincount`.
- Location/tenure redistribution maps are stored as dense index/weight tensors
  rather than Python dictionaries.
- Experimental Numba scatter kernels were added but are **off by default**
  (`P.use_numba_scatter=false`) because micro-kernel dispatch was slower in the
  full equilibrium run.

Warm one-iteration benchmark:

- Original one-iteration distribution timer: `4.214s`
- After dense/bincount distribution path: `3.110s`
- Moment output unchanged to numerical precision.

Accepted-equilibrium benchmark with the current path:

- file: `python_fast_ge_converged_distopt_dense.json`
- distribution timer: `103.30s` over `26` calls
- total wall time: `277.42s`
- moments still match MATLAB with max absolute difference `0.00531`

The full-run timing was unstable in this session and slower than the earlier
baseline even though the isolated one-step distribution timer improved. Treat
the current code as a correctness-preserving first pass, not the final
distribution solution. The next performance step should be a single compiled
one-period distribution transition kernel, not more small scatter wrappers.

## Distribution Optimization Pass 2

Changes:

- Added a compiled fast-stat transition kernel for the repeated GE distribution
  calls.
- The kernel propagates each mass cell directly through location choice, tenure
  transactions, savings interpolation, and child aging.
- The full-stat/event-study path remains on the transparent Python
  implementation, so final diagnostics are still computed by the original
  translated logic.

Warm one-iteration benchmark:

- Previous dense/bincount one-iteration distribution timer: `3.110s`
- Compiled fast transition distribution timer: `1.080s`
- Moment output is identical to the previous Python implementation up to
  `4.6e-13` in the comparison table.

Direct theta benchmark impact:

- Previous cached-eval direct-theta total wall time: `250.39s`
- Compiled fast transition direct-theta total wall time: `196.57s`
- Previous Python distribution timer: `94.48s` over `26` calls
- Current Python distribution timer: `30.14s` over `26` calls
- MATLAB distribution timer for the same direct-theta benchmark: `41.52s` over
  `26` calls

## Howard Eval Optimization Pass 1

Changes:

- Added two fused Numba kernels in `dt_cp_model/kernels.py`:
  `eval_renter_block_kernel` and `eval_owner_block_kernel`.
- The renter kernel collapses the eval-mode renter block (Cobb-Douglas
  consumption-housing aggregator with renter cap, plus continuation-value
  interpolation at `stored_bp`) into one tight loop over `(b, c)` cells with
  no intermediate `(Nb, nc)` allocations.
- The owner kernel does the equivalent collapse for each owner tenure level
  with deterministic housing services `chi * H_own[ten-1]`.
- The kernels are gated on `NUMBA_AVAILABLE` and a new `P.use_eval_kernel`
  switch that defaults to true; the original NumPy-vectorized eval block
  remains as a fallback for environments without Numba.

Warm direct-theta benchmark, eval-kernel toggle isolation in the same
process (so all earlier Numba kernels are already cached):

- `use_eval_kernel=False`: total `68.15s`, Howard eval `27.29s` over `18`
  calls.
- `use_eval_kernel=True`: total `51.80-58.67s` typical (run-to-run noise),
  Howard eval `21.36-23.94s` over `18` calls.

Updated direct-theta benchmark with the eval kernel on:

- file: `solve_theta_x0_fast_eval_kernel.json`
- Python total wall time: `~52-58s` (run-to-run, system noise).
- Python Howard eval time: `~22-25s` over `18` calls.
- Python full Bellman time: `~21-25s` over `8` calls.
- Python distribution time: `~8-11s` over `26` calls.
- MATLAB total wall time: `268.83s` (unchanged reference benchmark).
- Max absolute moment difference vs MATLAB: `0.00547`.
- Prices differ by less than `7e-5`.
- Distribution mass is conserved to floating-point precision.

Interpretation:

- The earlier `196.57s` Python total appears to have included Numba JIT
  compilation cost; warm-cache repeated runs at the same setup land closer to
  `60-70s` even before the eval kernel.
- Holding cache state fixed, the new eval kernel removes about `~5-6s` from
  Howard eval per direct-theta run and brings Python eval below MATLAB on a
  per-call basis when measured in the same warm session.
- Total Python direct-theta runtime is now `~5x` faster than the matching
  MATLAB benchmark with parity preserved at `0.00547`.
- Remaining eval-mode time is dominated by the post-eval tenure-choice,
  location-logit, and fertility-logit blocks, which still run in NumPy in
  both eval and full Bellman; that is the next candidate for compilation.

## Howard Eval Optimization Pass 2

Compiled the per-`j` post-Bellman blocks that previously ran in pure NumPy in
both eval-mode and full-mode:

- `tenure_choice_kernel` (in `dt_cp_model/kernels.py`): replaces the nested
  `(id_, to, tn, nn, cs, b)` loop and its `interp_on_grid` /
  `interp_vector` calls. Handles all three branches (move-while-renter, owner
  buys with possible birth-grant, sell-and-rent) in a single Numba loop,
  including `birth_dp`, `birth_entry_grant`, and `dp_arr` / `bmo`
  feasibility masks. Toggle: `P.use_tenure_kernel` (default true).
- `location_logit_kernel`: stable-logsumexp location choice across the
  `I=2` location set with `iidx`/`iwt` interpolation for the move case.
  Toggle: `P.use_loc_kernel` (default true).

Parity:

- vs prior Python output (iter6 warm): max abs diff `8.88e-16` — bit-exact
  up to floating-point rounding.
- vs MATLAB direct-theta x0 fast: max abs diff `0.00547` — unchanged.

Warm direct-theta benchmark with both new kernels enabled:

- file: `solve_theta_x0_fast_post_kernel.json`
- Python total wall time: `~34-37s` (run-to-run, system noise).
- Python Howard eval time: `~7.6-9.5s` over `18` calls (`0.42-0.53s` per call).
- Python full Bellman time: `~15-19s` over `8` calls (`1.9-2.4s` per call).
- Python distribution time: `~8-11s` over `26` calls (`0.31s` per call).
- MATLAB total wall time: `268.83s` (unchanged reference benchmark).
- Speedup vs MATLAB: `~7-8x`.

Cumulative effect of all warm-cache Bellman optimizations on the same
direct-theta benchmark, with the Numba binary cache held fixed:

| Configuration | Total | Per-eval | Per-full |
|---|---:|---:|---:|
| no Howard kernels (baseline) | `68.15s` | `1.52s` | `3.54s` |
| + eval renter/owner kernel | `~52-58s` | `~1.33s` | `~3.10s` |
| + tenure choice kernel | `~36-37s` | `~0.60s` | `~2.5s` |
| + location logit kernel | `~34-37s` | `~0.45s` | `~2.0s` |

Per-call breakdown at the final configuration:

- Per-eval call: `0.45s` (down from `1.52s` warm baseline, `~3.4x` faster).
- Per-full call: `1.96s` (down from `3.54s`, `~1.8x` faster).
- Per-dist call: `0.31s` (unchanged from the earlier distribution pass-2).

Remaining bottleneck:

- The fertility logit block at the bottom of each `j` step is still NumPy,
  but it is small (`~100us` per `j`).
- The full Bellman inner block (`solve_bellman_core`'s `not eval_mode`
  branch) still iterates the renter / owner golden-section per `(i, ten)`
  with NumPy housekeeping around the kernel. Full Bellman per call
  (`~2.0s`) now exceeds eval per call by `~4x` rather than `~2x` because
  eval is so much faster.
- `apply_child_aging` is one BLAS matmul per parity index; not a target.

## Howard Eval Optimization Pass 3 (2026-04-26)

Compiled the full-Bellman renter and owner inner blocks (the `not eval_mode`
branches) end-to-end:

- `full_renter_block_kernel`: replaces the `for c in range(nc)` loop +
  golden-section dispatch + post-search NumPy block (`surplus_nc`,
  `ct_nc`, `ht_nc`, cap mask, `co_nc`, `ho_nc`, infeasibility override)
  with one Numba kernel that loops over `(c, b)` and inlines both the
  golden-section search and the post-search consumption / housing
  bookkeeping. Toggle: `P.use_full_kernel` (default true).
- `full_owner_block_kernel`: same idea for each owner tenure `ten >= 1`.
  Each `(i, ten)` becomes one Numba call instead of `nc=24` Numba
  dispatches plus a NumPy post-block.

Parity:

- vs prior Python (iter6 warm): max abs diff `8.88e-16` — bit-exact.
- vs MATLAB direct-theta x0 fast: max abs diff `0.00547` — unchanged.

Warm direct-theta benchmark with all four passes enabled:

- file: `solve_theta_x0_fast_full_kernel.json`
- Python total wall time: `27.05s` (run-to-run noise gives `26-29s`).
- Python Howard eval time: `~7.0s` over `18` calls (`0.39s` per call).
- Python full Bellman time: `~10.7s` over `8` calls (`1.33s` per call).
- Python distribution time: `~7.0s` over `26` calls (`0.28s` per call).
- MATLAB total wall time: `268.83s` (unchanged reference benchmark).
- Speedup vs MATLAB: `~10x`.

Updated cumulative table:

| Configuration | Total | Per-eval | Per-full |
|---|---:|---:|---:|
| no Howard kernels (baseline) | `68.15s` | `1.52s` | `3.54s` |
| + eval renter/owner kernel | `~52-58s` | `~1.33s` | `~3.10s` |
| + tenure choice kernel | `~36-37s` | `~0.60s` | `~2.5s` |
| + location logit kernel | `~34-37s` | `~0.45s` | `~2.0s` |
| + full-mode renter/owner kernel | `~27s` | `~0.39s` | `~1.33s` |

What's left in `solve_bellman_core` after pass 3:

- The outer `for j` loop is sequential by construction (backward
  induction).
- The `for i` and `for ten` loops are `I=2` and `nt=4`; cheap.
- The fertility logit block (`~100us` per `j`).
- `apply_child_aging` (BLAS matmul, already fast).
- Per-call setup (`dp_arr`, `bmo`, `hcost`, `heq`, `Vbq`) runs once per
  call and is small.

## Real-Grid Benchmark (`benchmark` setup, Nb=80) — 2026-04-26

All speedup numbers above are at the small `fast` setup (`Nb=60`). The real
calibration grid uses `benchmark` (`Nb=80`). At the real grid the
speedup is much smaller because MATLAB's per-Bellman cost grows more slowly
with `Nb` than Python's — the gap that was 10x at `Nb=60` shrinks.

Direct-theta x0 at `benchmark` setup, `max_iter_eq=120`:

| | MATLAB | Python (cold) | Python (warm) |
|---|---:|---:|---:|
| `fast` (Nb=60) | `268.83s` | `41.1s` | `28.9s` |
| **`benchmark` (Nb=80)** | **`152.0s`** | **`81.7s`** | **`73.4s`** |

Apples-to-apples speedup at the real grid: **~2.1x** (`152.0s` MATLAB vs
`73.4s` Python warm), not ~10x. Both numbers are at the same
`max_iter_eq=120` and the same `theta=x0`.

The `fast` setup remains a useful unit-test harness for parity checks but
should not be cited as the calibration-grid speed.

Bridge-bench theta with slide-bench overrides (`E_C=0.1318`,
`r_bar_C=0.0762`, `H_own` from saved MAT, `hR_max=5.1`):

- Python at `benchmark` setup: `~90s` per solve, warm.
- Python total SMM loss against the live target table (no geography
  inversion): `30.80`. MATLAB `hard_loss` saved in
  `repair_candidate_..._localbaseline.mat` is `52.87`, but that includes
  the geography-inversion penalty, which we skip.
- All 17 hard-target moments match MATLAB within `±0.04` absolute.

## Parallel Kernels Pass (2026-04-26)

Added `@njit(cache=True, parallel=True)` and `for c in prange(nc):` to the
four block kernels:

- `eval_renter_block_kernel`
- `eval_owner_block_kernel`
- `full_renter_block_kernel`
- `full_owner_block_kernel`

Each `c` iteration writes to disjoint output cells, so the parallelization
is race-free. The other kernels (`tenure_choice_kernel`,
`location_logit_kernel`, `forward_distribution_fast_kernel`) remain
single-threaded for now.

Parity:

- vs prior Python (iter6 warm): max abs diff `8.88e-16` — bit-exact.
- vs MATLAB direct-theta x0 at `benchmark` setup: max abs diff `0.00282`.

Direct-theta x0 at `benchmark` setup, warm:

| | MATLAB | Python (serial) | Python (parallel) |
|---|---:|---:|---:|
| Total | `152.0s` | `73.4s` | `54.4s` |
| bellman_full per call | `9.76s` | — | `1.14s` |
| bellman_eval per call | `0.79s` | — | `0.97s` |
| distribution per call | `0.67s` | — | `0.47s` |
| Speedup vs MATLAB | — | `2.07x` | **`2.79x`** |

CPU utilization during the run averaged `165-167%` of wall (i.e. `~1.65`
effective cores busy on an 8-core machine), not the ideal `8x`. The
`prange(nc=24)` is the bottleneck: small loop body in eval-mode kernels
plus thread-spawn overhead. The serial parts (`j` loop, GE loop, forward
distribution, tenure choice, location logit) are still single-threaded.

Best file: `benchmarks/solve_theta_x0_bench_parallel.json`.
