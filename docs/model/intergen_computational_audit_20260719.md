# Computational Audit of the Intergenerational Housing-Fertility Model

**Date:** 2026-07-19

**Production object audited:** `code/model/intergen_housing_fertility/` under the M5 contract

**Production code modified:** none

**Audit artifacts:** `output/model/computational_audit_20260719/`

## Executive verdict

The default M5 numerical path is computationally coherent and well protected by
invariants, but it is slower and harder to maintain than necessary. There is no
evidence that a wholesale C++ or Rust rewrite is warranted. The active Bellman
blocks are already compiled with Numba. The largest remaining costs are the
Python-level Markov forward-distribution path and the number of times the full
fixed-price household problem is solved while finding the market-clearing
price.

The current strict production record takes `47.44` seconds per M5 evaluation on
Torch. Of that:

| Stage | Seconds | Share |
|---|---:|---:|
| Damped price iterations | 31.51 | 66.4% |
| Scalar Brent refinement | 12.83 | 27.0% |
| Final full-statistics solve | 3.09 | 6.5% |
| **Total** | **47.44** | **100%** |

Thus 93.5% of the production solve is price-search work. The local function
profile, though run under heavy host contention, gives a stable decomposition
of that work: across 34 fixed-price solves, the Markov forward pass accounts
for 64.3% of cumulative time and the Bellman pass for 35.6%.

The complete canonical 15-row target/model/gap/weight/contribution table remains
`output/model/intergen_income_disciplined_recalibration_20260716/report/target_fit_full.csv`;
this computational audit does not replace or reinterpret that calibration
report.

The audit also found two confirmed correctness defects in supported but
nondefault paths and one nondefault measurement defect:

1. changing income persistence or same-sized income weights through the generic
   override API can leave the transition matrix stale;
2. the Python owner fallback decodes the flattened family/child-state index in
   the wrong order, affecting the documented monotone-cubic path;
3. the old-age nonhousing median uses a collapsed-income denominator when
   retirement income is allowed to vary with the persistent income state.

None of these three defects is active in the default linear-interpolation M5
solve. They should nevertheless be fixed before the affected switches are used
for quantitative work.

## Scope and method

The canonical status file identifies M5 as the current provisional paper
calibration. E1 is an experimental fork, so this audit treats
`intergen_housing_fertility` as production and discusses the forks only as a
maintenance issue.

The audit combined:

- an exact fixed-price M5 reproduction at the collected equilibrium price;
- a full tight GE profile at `J=17`, `Nb=120`, five income states, five owner
  rungs, `max_iter_eq=40`, and `tol_eq=2.5e-5`;
- separate fast-statistics and full-statistics fixed-price profiles;
- one- versus four-thread and Python- versus compiled-scatter probes;
- static review of the solver, kernels, parameters, calibration, optimizer, and
  tests;
- targeted reproductions of suspected correctness failures;
- the complete production-package test suite.

Local profiles used Python 3.10.10, NumPy 1.24.3, Numba 0.58.1, and one Numba,
OMP, MKL, and OpenBLAS thread unless otherwise stated. The laptop carried a
load average between roughly 4 and 15 during the runs. Absolute local wall times
are therefore not portable; function shares, call counts, exact-output checks,
and the canonical Torch timing are the relevant evidence.

## 1. Profiling results

### 1.1 Full strict GE solve

The local tight profile made 34 calls to `solve_markov_income_at_prices`:

- 18 damped fixed-point iterations;
- 15 scalar-refinement evaluations;
- 1 final full-statistics solve.

The local equilibrium is the already documented macOS platform reproduction:
loss `8.9989694`, residual `9.32e-6`, versus the canonical Torch loss
`9.0444221`, residual `1.48e-5`. The local and Torch prices differ by only
`2.8e-5`, but kinked/discrete decisions move the most sensitive moment by
`0.00952`. Both equilibria satisfy the declared tolerance. This is a
cross-platform equilibrium-tolerance issue, not run-to-run local noise; the
saved-price mapping reproduces the Torch moments exactly.

The full local profile decomposes as follows:

| Function/block | Cumulative seconds | Share of local GE profile |
|---|---:|---:|
| `forward_distribution_markov_income` | 124.91 | 64.3% |
| `solve_bellman_full_markov_income` | 69.10 | 35.6% |
| `refine_one_market_markov_income` | 59.41 | 30.6% |
| `realize_current_cross_section` | 53.79 | 27.7% |
| Python `scatter_redistribute` calls | 53.61 | 27.6% |
| Compiled owner blocks | 36.28 | 18.7% |
| Compiled tenure logits | 12.44 | 6.4% |
| Compiled renter blocks | 9.68 | 5.0% |

The entries overlap: the realization and scatter rows are components of the
forward pass, and owner/renter/tenure rows are components of the Bellman pass.

### 1.2 One fast fixed-price evaluation

The `fast_stats=True` profile approximates the object repeated during market
clearing. At one thread and with the production Python scatter path:

| Block | Seconds | Share |
|---|---:|---:|
| Markov forward distribution | 6.06 | 60.5% |
| Bellman solve | 3.93 | 39.3% |
| Total | 10.01 | 100% |

The level is inflated by host contention, but the decomposition agrees with the
full GE profile.

### 1.3 Final full-statistics solve

The full fixed-price profile shows that reporting is appreciably more expensive
than market-clearing statistics:

| Block | Seconds | Share |
|---|---:|---:|
| Markov forward distribution plus reporting setup | 20.01 | 81.5% |
| Bellman solve | 4.52 | 18.4% |
| `assign_current_cross_section_to_beginning_assets` | 8.86 | 36.1% |
| `compute_markov_statistics` | 0.98 | 4.0% |

The beginning-asset/current-choice reassignment, not the scalar moment formulas,
is the expensive reporting object. This is paid only once per equilibrium, so
it is secondary for calibration throughput but important for interactive
diagnostic generation.

### 1.4 Existing compiled scatter and threads

The existing `use_numba_scatter=True` switch reproduced all 15 saved-price
moments to `1.7e-14`. Under the contended local benchmark it reduced a
one-thread fast evaluation from `10.01` to `8.82` seconds, about 12%. Raising
Numba from one to four threads with compiled scatter reduced it further to
`6.82` seconds, about 23% relative to the one-thread compiled-scatter case.

These are promising probes, not production acceptance tests. The cluster runs
one-core candidate evaluations and already obtains parallelism across chains and
candidates; using four threads per solve could reduce total throughput if CPU
allocations are unchanged. The compiled-scatter switch should be tested over a
representative GE parameter panel before being enabled by default.

## 2. Correctness review

### 2.1 What passed

At the canonical saved M5 price:

- all 15 model moments match the collected Torch table to
  `1.8e-15` with the production scatter path;
- total distribution mass is `0.99999999999995`;
- the distribution has no negative or nonfinite mass;
- tenure, location, and fertility probabilities remain in `[0,1]` and add to
  one on active rows, with maximum tenure float32 sum error `6.4e-8`;
- savings policies are finite and no policy exceeds the wealth grid by more
  than the golden-section tolerance;
- raw value-function monotonicity violations occur on about `0.01%` of adjacent
  alive-grid pairs, but only `8.8e-11` of distribution mass is adjacent to
  them; they are concentrated near penalty/dead-state transitions;
- market residual is `1.48e-5`, below the declared `2.5e-5` tolerance.

The existing code has valuable defenses: explicit dead-state mass gates and
censuses, mass-conservation assertions, current-choice versus beginning-asset
timing separation, strict candidate-selection gates, exact contract metadata,
and focused debt, survival, bequest, transfer, and transaction tests.

### 2.2 Confirmed findings

#### C1. Stale income transition matrix after generic overrides

**Severity:** high for the generic API; inactive in M5's explicit production
helper.

`apply_overrides` rebuilds `Pi_z` only when it is absent or has the wrong shape.
If `income_shock_persistence` changes while the number of income states stays
fixed, the old shape-compatible matrix is retained and merely row-normalized.
Changing same-length `z_weights` has the same problem.

Targeted reproduction:

- changing persistence from `0.85` to `0.10` left `Pi_z` bit-identical;
- changing weights to `[0.8,0.1,0.1]` left a matrix whose invariant-weight
  error is `0.075`.

Production M5 is protected because `income_process_overrides` passes `Pi_z`
explicitly. Direct comparative statics through the public override API are not.

**Recommendation:** if persistence or weights change and `Pi_z` is not
explicitly supplied in the same override, rebuild it. If `Pi_z` is explicit,
validate row sums and either validate or derive its stationary distribution.

#### C2. Wrong flattened-state decoding in the Python owner fallback

**Severity:** high for the fallback; inactive in the default compiled-linear
profile.

`flat_nc` uses Fortran order,

\[
c=n+n_{parity}cs,
\]

so parity varies fastest. The Python owner fallback instead divides by the
number of child states. With the production `3 x 4` state layout, 10 of 12
columns select the wrong `(n,cs)` collateral floor. The error appears in both
the Markov and inherited non-Markov fallback blocks.

The default linear M5 solve uses `full_owner_block_kernel` and is unaffected.
The documented `interp_method="monotone_cubic"` Markov path deliberately falls
back to Python, so that sensitivity path is unsafe whenever `phi` or another
financing hook differs across family/child states.

**Recommendation:** decode with `n = c % n_parity` and
`cs = c // n_parity`, then add compiled-versus-fallback array parity tests under
state-dependent `phi`.

#### C3. Collapsed-income denominator in a nonlinear old-age median

**Severity:** medium; inactive when `retirement_income_z_scale=0`, as in M5.

`compute_markov_statistics` first collapses the income-state dimension and
then calls `compute_statistics`. Several nonlinear renter and wealth objects are
recomputed correctly on the full income-resolved distribution, but
`old_nonhousing_wealth_to_income_median_6575` remains calculated with income
evaluated at `z=1`.

**Recommendation:** recompute this median and all related nonlinear
wealth-to-income quantiles directly over `(b,z)` cells whenever retirement
income depends on `z`. Add a test with nonzero retirement-income loading.

#### C4. Production target system is assembled in two places

**Severity:** medium source-of-truth risk.

The calibration registry named
`candidate_replacement_income_disciplined_v1` contains 14 moments. The actual
M5 runner wraps it with `target_system` and appends
`aggregate_mean_occupied_rooms_18_85`, producing the audited 15-moment system.
Some report drivers remember to append the same row; generic callers of
`get_target_set` do not.

This split initially caused the audit harness to compute loss `8.9753` instead
of the correct saved-price loss `9.0444`, even though the model moments were
identical. The production M5 chain is correct, but the API makes accidental
underidentification or understated loss too easy.

**Recommendation:** put the complete 15-moment M5 target system in one
canonical registry and make every runner consume it without mutation.

#### C5. Weak parameter-domain and override validation

**Severity:** medium robustness risk.

The dynamic `SimpleNamespace` override path accepts unknown keys and many
invalid domains. Direct callers can supply misspelled parameters, invalid
preference curvature, nonpositive logit scales, invalid housing menus, or
degenerate grid boundaries without an early, informative failure. Search bounds
protect current calibration jobs, but the public solver boundary does not.

**Recommendation:** add a centralized `validate_parameters` gate immediately
after overrides. Reject unknown keys except an explicit set of experiment hooks.
This is more valuable than converting the entire project to a new type system.

### 2.3 Reproducibility qualification

The macOS/Torch equilibrium difference is already documented in
`output/model/intergen_m5_draft_refresh_20260717/gate0/PLATFORM_NOTE.md`.
The current audit reproduces it exactly. The household and KFE mapping at a
fixed price is cross-platform stable to reported precision; the small difference
comes from the equilibrium search accepting different nearby points under a
kinked demand schedule and tolerance `2.5e-5`.

For levels, the canonical Torch report remains the source of truth. For policy
deltas, both benchmark and counterfactual must be solved in the same process and
environment. If cross-platform levels must agree more tightly, the price root
contract needs a tighter and more reproducible stopping rule; a language rewrite
would not by itself solve this issue.

## 3. Test review

| Check | Result |
|---|---|
| Python compile check for the audit driver | pass |
| `unittest` production suite | 64/64 pass in 3.95 seconds |
| Top-level pytest-style functions, invoked directly | 9/9 pass |
| Total test functions exercised | 73 pass |
| Normal `pytest` collection | fails with exit 139 and no Python traceback |
| Exact fixed-price M5 reproduction | pass, all 15 moments within `1.8e-15` |
| Compiled-scatter saved-price equivalence | pass, all 15 moments within `1.7e-14` |

The pytest crash is an environment/test-runner defect even though the same test
bodies pass when invoked through `unittest` or directly. It prevents dependable
CI and should be diagnosed separately. Until then, the repository should expose
one stable command that runs both the 64 class-based and 9 top-level tests.

Highest-value missing tests:

1. persistence/weight override consistency and invariant-distribution checks;
2. compiled versus Python Bellman parity with state-dependent financing terms;
3. full income-resolved old-age nonlinear moments under retirement-income
   dispersion;
4. one named M5 end-to-end golden regression covering moments, residual, mass,
   probability adding-up, and boundary mass;
5. compiled-scatter equivalence across several GE parameter points;
6. a cross-platform price-root tolerance test based on demand/supply residuals,
   not bitwise equality.

## 4. Code clarity, length, and style

The production package contains 12,276 lines of Python excluding tests and
14,466 including its tests. `solver.py` alone is 6,287 lines with 102 top-level
functions. The largest functions are:

| Function | Lines | Approximate branch nodes |
|---|---:|---:|
| `compute_statistics` | 581 | 125 |
| `solve_bellman_core` | 508 | 85 |
| `forward_distribution` | 484 | 84 |
| `forward_distribution_markov_income` | 462 | 84 |
| `solve_equilibrium` | 398 | 76 |
| `solve_bellman_full_markov_income` | 360 | 58 |

The main clarity problems are structural rather than cosmetic:

- Markov and non-Markov Bellman/forward/statistics paths duplicate indexing and
  accounting logic. The fallback-index and collapsed-moment defects are examples
  of the drift this permits.
- `SimpleNamespace` plus unrestricted overrides obscures the model contract and
  makes typos silently legal.
- `calibration.py` combines target values, target-object documentation,
  diagnostic search routines, and result extraction.
- `local_panel.py` combines candidate generation, optimization algorithms,
  persistence, acceptance, and reporting.
- numerical conventions such as dead-value sentinels, interpolation clipping,
  golden tolerance, axis order, and small-mass thresholds are repeated across
  Python and Numba implementations.
- `production_profile.py` still names the July 9 repair target system, and
  `IMPLEMENTATION_STATUS.md` still says formal calibration is not implemented,
  despite M5 being live.

Experimental isolation has also produced three near-complete package copies:

| Package | Source lines |
|---|---:|
| Production intergen | 12,276 |
| Sequential-fertility fork | 12,535 |
| Equivalence-scale/sequential fork | 13,153 |

Keeping E1 isolated while it is experimental was defensible. If it survives,
the shared numerical core should be factored once rather than maintaining about
38,000 copied source lines. This should happen after the experiment's state
architecture stabilizes, not during the current calibration cycle.

Positive style features include explicit model-object comments, named
production configurations, deterministic seeds, clear checkpoint records,
well-documented target objects, and comments that distinguish intended
shortcuts from model claims. The main numerical functions are readable locally;
their excessive scope makes global reasoning and parity maintenance difficult.

## 5. Prioritized speed recommendations

### P1. Reduce the number of fixed-price solves

This has the largest ceiling because price work consumes 93.5% of the canonical
runtime. The present one-market solver first performs damped fixed-point
iterations and then starts a separate bracketed scalar refinement. The audit's
local solve used 33 fast price evaluations before the final report solve; the
canonical Torch solve used 28.

Recommended experiment:

1. construct a deterministic price bracket directly;
2. solve excess demand with safeguarded Brent/Illinois from the bracket;
3. cache already evaluated prices and their solutions;
4. use the current damped loop only as a fallback when no bracket is found.

If the fast evaluation count fell from 28 to 12 with unchanged per-evaluation
cost, the canonical solve would fall from about 47 seconds to about 22 seconds,
roughly `2.1x`. This is an arithmetic scenario, not yet a benchmark. The new
root solver must be tested over a representative calibration panel because
housing demand is kinked and need not be globally monotone.

### P2. Fuse and compile the Markov forward path

The Markov forward pass is 64% of tight-GE runtime. It still loops in Python over
age, income, tenure, family size, and child state, with millions of small
scatter/reduction operations and repeated temporary arrays. The inherited
non-Markov path already demonstrates a fused compiled forward kernel.

A fused Markov kernel should combine current-choice realization, tenure
transactions, savings interpolation, income transitions, child aging, and mass
accounting in one or a few compiled passes. If the forward pass were made `2x`
faster, Amdahl's law implies about `1.47x` end-to-end speed; at `3x`, about
`1.75x`, holding price-evaluation count fixed.

The existing compiled-scatter switch is a useful interim step. Its measured
gain was only about 12% for the fast solve because Python loop and allocation
overhead remains.

### P3. Reuse the accepted-price household solution

The price solver evaluates the accepted or nearly accepted price with
`fast_stats=True`, then recomputes the Bellman solution from scratch for the
full report. Retaining the value/policy arrays for the best price and recomputing
only the full distribution/reporting objects would remove one Bellman pass.
The ceiling is modest—roughly 2–4% of total canonical time—but the change is
conceptually contained.

### P4. Add carefully verified policy warm starts

The Markov path has no Howard/policy-evaluation reuse. Adjacent prices and
adjacent optimizer candidates often have similar policies. A previous-price
savings policy could narrow or seed golden searches, but must retain an
unrestricted fallback and a dense/globality audit. Past diagnostics found
nontrivial branch-level globality gaps even when mass-weighted gaps were small,
so an unverified local bracket would trade correctness for speed.

### P5. Optimize full reporting after the GE loop

For interactive packets, the current-choice/beginning-asset reassignment is the
largest final-reporting cost. It can be fused with the realized-distribution
pass or compiled independently. This will improve diagnostics and
counterfactual packet generation more than calibration throughput.

### P6. Use threads only where the resource contract supports them

Four Numba threads improved one local fast solve by about 23% relative to the
one-thread compiled-scatter run. The current Torch strategy gets parallelism
from independent chains/candidates and generally allocates one CPU per task.
Threading should be benchmarked as throughput per allocated core, not latency of
one solve, before changing cluster defaults.

## 6. Implementation order

Recommended order if the project chooses to act:

1. fix C1–C4 and add their regression tests;
2. provide one stable test command and resolve the pytest crash;
3. run a representative GE A/B panel for `use_numba_scatter=True`;
4. prototype the direct bracketed price solver in a nonproduction module and
   compare price, residual, all moments, and solve counts over many candidates;
5. prototype a fused Markov forward kernel with exact mass and moment parity;
6. only then consider policy warm starts or structural code splitting;
7. factor shared production/E1 numerical code after E1's architecture is
   accepted.

The likely combined gain from fewer price evaluations plus a fused forward pass
is several-fold, while retaining Python for orchestration and reporting. That is
the appropriate optimization frontier before considering a C++ or Rust port.

## Reproduction commands

From `code/model`:

```bash
PYTHONPATH=$PWD NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python tools/profile_intergen_computation.py --mode fixed-price-fast

PYTHONPATH=$PWD NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python tools/profile_intergen_computation.py --mode fixed-price

PYTHONPATH=$PWD NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python tools/profile_intergen_computation.py --mode tight-ge
```

The raw `.pstats`, sorted text profiles, summaries, and invariant packets are in
`output/model/computational_audit_20260719/`.
