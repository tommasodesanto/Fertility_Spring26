# Optimized Intergenerational Housing-Fertility Model

> **Experimental physical refactor.** This package is isolated from
> `intergen_housing_fertility` and no production runner imports it. The source
> package remains unchanged. Promotion requires the parity and strict-repeat
> gates below plus a representative parameter-panel benchmark.

This package is a correctness-first refactor of the production M5 model. Its
current accepted changes are:

- an atomic, immutable 15-moment M5 target system;
- strict override validation and consistent income-process rebuilding;
- centralized Fortran-order family-state layout helpers;
- income-resolved retirement wealth/income statistics;
- compiled Markov scatter operations enabled by default;
- a safeguarded directional one-market Brent search with legacy fallback;
- price-evaluation caching and reuse of the accepted Bellman solution;
- a package-owned exact M5 configuration and reproducible parity benchmark.

An attempted fully fused Markov age kernel was numerically accurate but slower
than the accepted compiled-scatter path, so it was measured, documented, and
removed rather than left as a dormant branch.

Run the complete stable test set from `code/model`:

```bash
PYTHONPATH=$PWD .venv/bin/python -m intergen_housing_fertility_optimized.tests.run_all
```

Run exact saved-price production parity:

```bash
NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
PYTHONPATH=$PWD .venv/bin/python -m intergen_housing_fertility_optimized.benchmark_m5 \
  --package both --mode fixed-price
```

Run the six-case fixed-price perturbation panel:

```bash
NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
PYTHONPATH=$PWD .venv/bin/python -m intergen_housing_fertility_optimized.parity_panel \
  --output ../../output/model/intergen_optimized_refactor_20260719/fixed_price_panel.json
```

Run the isolated strict optimized equilibrium:

```bash
NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
PYTHONPATH=$PWD .venv/bin/python -m intergen_housing_fertility_optimized.benchmark_m5 \
  --package optimized --mode tight-ge
```

See `REFACTOR_REPORT.md` for current evidence, limitations, and the promotion
checklist.

The complete local/Torch promotion battery is specified in
`PROMOTION_RUNBOOK.md`. It covers every live M5 parameter bound, the configured
price bracket, actual and forced-fallback root cases, strict repeats, all
standard diagnostic packet inputs, and a checkpointed ten-candidate
calibration-throughput smoke. The Torch collector now passes all gates. The
package remains non-production and no existing import has been redirected,
because numerical promotion and production adoption are separate decisions.

The bounded six-hour continuation workflow is documented in
`CALIBRATION_RUNBOOK.md`. It preserves the M5 specification and target system,
uses recoverable one-CPU Torch chains, and promotes only twice-repeated strict
winners; loose search losses are never compared with canonical M5.

The invalid July 22 draft ledger is retained in the status history, but the
current `new_moment_profile.py` no longer equates CEX estimates with structural
parameters. It reproduces four model-feasible CEX auxiliary regressions and a
matched PSID four-year tenure Brier score on simulated households. A local
29-solve transformed-coordinate audit is finite and full rank 14/14
(condition number 858); the packet is under
`output/model/intergen_new_moment_jacobian_20260722/full/`. The three-hour Torch launchers are
`cluster/submit_three_hour_new_moments.sh` and
`cluster/submit_three_hour_new_moment_collector.sh`. This is a provisional
one-shot calibration profile, not a replacement for M5 until a strict winner
is collected and reviewed.

## Inherited model documentation

This package is the new one-market quantitative implementation for the
intergenerational housing allocation and fertility project. It starts from the
active workhorse lifecycle code under `code/model/dt_cp_model/`, but keeps that
code path unchanged.

For the simplest current-model run, open and run:

```bash
code/model/run_intergen_model.py
```

That file lists the current source theta, target set, grid, and output folder
at the top.

In Spyder, open `code/model/run_intergen_model.py` and press Run. The file runs
the model through `code/model/.venv/bin/python` without replacing the Spyder
kernel, then loads `solution_summary`, `moments`, `target_fit`, `age_profiles`,
`room_bin_fit`, `first_look_density_path`, `first_look_total_wealth_path`,
`first_look_total_wealth_density_path`, `solution_cache_path`,
`first_look_policy_lines`, `first_look_market_summary`, and the total-wealth
density/support tables for inspection in the Variable Explorer.
The runner also prints macOS `CPU_Speed_Limit` before and after the solve. If a
slow run has a low speed limit, the same model path is being throttled by the
laptop rather than taking a different equilibrium path.
After one successful solve, set `REFRESH_PLOTS_FROM_SAVED_SOLUTION = True` in
the runner to rebuild the graph packet from `solution_cache.pkl` without
solving the model again. Leave `FAST_REFRESH_FROM_SAVED_SOLUTION = True` for
quick first-look plot iteration; it skips the larger standard-diagnostic,
contact-sheet, age-profile, room-bin, rung, and threshold-plot rebuilds while
still refreshing the first-look plots and CSVs.

Run from `code/model` when using the package CLI directly:

```bash
.venv/bin/python -m intergen_housing_fertility.cli smoke --quiet
.venv/bin/python -m intergen_housing_fertility.cli solve --max-iter-eq 20 --quiet
.venv/bin/python -m intergen_housing_fertility.cli diagnostics --fixed-prices --outdir ../../output/model/intergen_housing_fertility_smoke_fixed --quiet
.venv/bin/python tools/build_intergen_mechanics_packet.py --max-iter-eq 10
```

The current pass uses one aggregate housing-services market, a 4-year decision
period, one dependent-child stage, a persistent Markov income state, lifecycle
income by age, owner housing rungs, continuous renter housing, down-payment
constraints, and transaction/sale wedges. The live diagnostic convention caps
renter housing at `hR_max=6.0` while the owner ladder remains
`H_own=[2,4,6,8,10]`; this keeps strong tenure/product segmentation while
allowing a thin first large-rental margin and reserving the upper 8/10-room
rungs for owners. A payment-to-income screen exists as an optional diagnostic
switch, but the default model leaves it off and uses the collateral/down-payment
constraint as the active finance restriction. It is not calibrated.

`IMPLEMENTATION_STATUS.md` is the live implementation record. Any future
simplification or deferred object should be added there in the same edit as the
code change.

`tools/build_intergen_mechanics_packet.py` is a non-production inspection
driver for the June 2026 one-market strand. It reads a saved theta, re-solves
the active intergen model, and writes standard diagnostics, target-fit tables,
room-bin/rung shares, simple and full first-look policies/market panels, the
same first-look policy panels against total wealth after the chosen tenure
branch, standalone aggregate liquid-wealth and total-wealth density plots with
childless/parent and renter/owner splits, wealth-support coverage diagnostics,
owner-rung diagnostics, owner-entry thresholds, a local solved-object cache, and
optional policy
proof-of-concept cases under
`output/model/intergen_mechanics_packet_YYYYMMDD/`.
