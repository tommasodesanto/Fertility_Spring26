# Python Port Plan For The April 2026 Discrete-Time Model

Status date: 2026-05-07

This subfolder is now the active implementation of the April 2026
center-periphery fertility model. The legacy MATLAB files used for parity
checks have been archived under
`calibration_archive/model_history_2026-05-07/legacy_matlab_2026-05-07/`.

Implementation status:

- Core setup, solver, distribution, statistics, and SMM objective are ported.
- The canonical Python model boundary is now `solve_theta(theta, setup_mode=...)`;
  SMM uses the same `apply_theta` path rather than reinterpreting parameter
  vectors separately.
- Optional Numba kernels accelerate the golden-section savings search when
  `numba==0.58.1` is installed.
- Reduced fast-setup parity benchmark passes with max absolute moment
  difference `0.00263` against MATLAB for a one-iteration GE check.
- Reduced fast-setup GE benchmark also matches MATLAB under the live solver
  stopping rule, with max absolute moment difference `0.00531` and prices
  within `5e-6`.
- Distribution optimization pass 2 compiles the fast-stat one-period mass
  transition. The isolated one-step distribution timer is now `1.08s` versus
  `4.214s` in the original Python path, with identical Python moments.
- Direct `solve-theta --setup fast --theta x0 --max-iter-eq 120` is the
  canonical parameter-vector model path. With cached Howard-eval interpolation
  and compiled fast distribution, it takes `196.57s` versus `268.83s` for the
  matching MATLAB direct-theta benchmark, with max absolute moment difference
  `0.00547`.
- Howard eval optimization pass 1 (2026-04-25) adds compiled
  `eval_renter_block_kernel` and `eval_owner_block_kernel` that fuse the
  eval-mode renter and owner blocks into one tight loop per
  `(i, ten)`. With the same Numba caches warm in a single session, this drops
  total direct-theta runtime to `~52-58s` and Howard eval to `~22-25s` over
  `18` calls, while preserving the `0.00547` max abs moment difference vs
  MATLAB. The earlier `196.57s` headline number was inflated by Numba JIT
  compilation in the same run; the apples-to-apples warm baseline is `~68s`.
- Howard eval optimization pass 2 (2026-04-25) adds compiled
  `tenure_choice_kernel` and `location_logit_kernel` for the per-`j` post-
  Bellman blocks that run in both eval and full mode. Warm direct-theta
  total drops to `~34-37s`, Howard eval to `~7.6-9.5s` (`0.45s` per call,
  down from `1.52s`), full Bellman to `~15-19s` (`1.96s` per call, down from
  `3.54s`). Speedup vs MATLAB direct-theta is now `~7-8x`. Parity vs MATLAB
  unchanged at `0.00547`. Both kernels are bit-exact vs prior Python at
  `8.88e-16`.
- Howard eval optimization pass 3 (2026-04-26) adds compiled
  `full_renter_block_kernel` and `full_owner_block_kernel` for the full-
  Bellman inner blocks. They eliminate the per-`c` Numba dispatch and the
  NumPy post-search arithmetic (consumption / housing / cap mask). Warm
  direct-theta total drops to `~27s`. Per-eval call is `0.39s`, per-full
  call is `1.33s`, per-dist call is `0.28s`. Speedup vs MATLAB direct-theta
  is now `~10x`. Parity vs MATLAB unchanged at `0.00547`. Bit-exact vs
  prior Python at `8.88e-16`.
- Full production-scale calibration has not been run from Python yet.

## Scope

Port the active discrete-time center-periphery lifecycle model:

- setup and calibration target construction
- fixed-price partial equilibrium and full housing-supply equilibrium
- backward Bellman recursion with tenure, location, and fertility logits
- forward distribution propagation
- live calibration moments and SMM loss
- benchmark/parity tools

Do not port archived search drivers, plotting batteries, or slide-specific
exporters in this first pass. They can be rebuilt against the Python solution
object after the solver is validated.

## Design Principles

1. Preserve economics first.
   The first Python solver uses the same state space and ordering conventions
   as MATLAB, including the one-shot completed-fertility architecture and
   financed-share interpretation of `phi`.

2. Make performance visible.
   Every benchmark records Bellman, distribution, and total wall-clock time.
   Speed changes are accepted only with moment parity checks.

3. Optimize high-return kernels only.
   The expensive path is Bellman plus distribution propagation. Diagnostics,
   plotting, and text exporters should not drive the translation effort.

4. Keep data layout explicit.
   Arrays use MATLAB-compatible dimensions:
   `(wealth, tenure, location, age, parity, child_state)`.
   Flattening over `(parity, child_state)` uses Fortran order where the MATLAB
   code does.

5. Separate parity from improvements.
   The initial solver mirrors the MATLAB golden-section architecture. Later
   acceleration can add Numba/JAX kernels, endogenous-grid variants, or
   batched parameter evaluation after parity is established.

## Implementation Steps

### Step 1: Package Skeleton

- Add `dt_cp_model` Python package.
- Require only NumPy initially, because the local default Python lacks SciPy and
  Numba. Optional acceleration dependencies can be added later.
- Provide a small CLI so benchmark commands are reproducible.

### Step 2: Setup Layer

- Translate `setup_parameters`.
- Translate `build_calibration_setup`.
- Support `benchmark` and `fast` setup modes.
- Implement MATLAB-style override semantics for:
  `kappa_fert`, `kappa_loc`, `eps_fert`, `eps_loc`, `H_own`, `H0`,
  `entry_shares`, income/pension fields, and grid dimensions.
- Rebuild child-transition matrices whenever fertility-state dimensions change.

### Step 3: Numerical Utilities

- Wealth grid construction.
- Linear interpolation with clipped extrapolation.
- Linear redistribution matrices represented as `(idx, weight)` pairs instead
  of sparse matrices.
- Stable log-sum-exp for location and fertility logits.
- Weighted medians/quantiles for room moments.

### Step 4: Bellman Solver

- Port the backward induction exactly:
  savings choice, tenure choice, location logit, fertility logit.
- Use vectorized golden-section search over wealth states for each
  `(parity, child_state)` column.
- Keep Howard policy evaluation:
  full solve for early/stalled iterations, cheap revaluation otherwise.
- Preserve key model details:
  renter cap at `hR_max`, owner housing services `chi * H_own`, birth
  down-payment grant logic, stochastic child aging, bequest utility, and
  additive DUE location utility.

### Step 5: Forward Distribution

- Replace MATLAB sparse redistribution matrices with cached interpolation maps.
- Use mass-preserving `numpy.add.at` redistribution for:
  location moves, tenure changes, savings, and child aging.
- Compute both fast equilibrium statistics and the full calibration moments.
- Preserve event-study housing moments:
  `housing_increment_0to1_eventstudy_t3` and
  `housing_increment_1to2_proxy_t3`.

### Step 6: Objective And Calibration

- Port the SMM objective.
- Preserve the current Stage A 13-parameter vector and legacy 12-parameter
  compatibility.
- Preserve the geography inversion over `E_C` and `r_bar_C`.
- Add a minimal random/PSO-ready objective wrapper only after a single solve is
  validated.

### Step 7: Benchmarks

Run three benchmark classes:

- Python smoke: tiny fixed-price solve, verifies arrays and finite moments.
- Python fast setup: small live model solve with timings.
- MATLAB parity: run the same reduced override in MATLAB and compare prices,
  TFR, ownership, population shares, migration, and housing moments.

The first accepted benchmark should report:

- wall time
- Bellman full/eval/distribution timing
- key moments
- maximum absolute and relative differences versus MATLAB where available

## Optimization Roadmap After Parity

1. Distribution kernel:
   replace repeated `numpy.add.at` calls with a Numba kernel or a custom Cython
   scatter-add.

2. Bellman full solve:
   move golden-section kernels to Numba, keeping NumPy interpolation parity.

3. Batched parameter evaluation:
   cache all price-independent tensors and evaluate candidate parameter points
   in a process pool.

4. Alternative savings solver:
   test endogenous-grid or monotone-policy methods only after the golden-section
   solver has passed parity. This model has tenure thresholds and borrowing
   constraints, so monotonicity must be verified before replacing the solver.

5. Objective search:
   use checkpoint-aware Python search with explicit run budgets, heartbeat
   summaries, and best-so-far artifacts, following the project long-run search
   safety rules.

## Validation Gates

No Python calibration run should be trusted until these pass:

- child-state transition rows sum to one
- probabilities are in `[0,1]` and sum to one over choices
- distribution mass is conserved before final normalization
- housing demand is finite and positive by location
- value/policy arrays are finite on reachable states
- Python and MATLAB reduced solves agree to tolerances documented in the
  benchmark output
