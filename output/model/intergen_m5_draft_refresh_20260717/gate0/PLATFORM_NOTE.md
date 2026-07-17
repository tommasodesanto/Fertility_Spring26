# Platform note: local M5 baseline vs. the canonical Torch record (2026-07-17)

## Canonical numbers

The **Torch record is canonical** for all reported M5 numbers: loss
`9.044422069071352`, residual `1.48e-05`, and the 15 moments in
`output/model/intergen_income_disciplined_recalibration_20260716/report/target_fit_full.csv`.
Nothing in this folder supersedes that record.

## What this folder contains

`m5_solution_cache.pkl` is the **local M5 baseline**: a fresh solve of the
exact M5 contract (theta from `results.json` winners.M5, arm `M5` via
`arm_contract`/`common_overrides` from
`code/model/tools/run_intergen_bequest_exit_chain.py`, `J=17`, `Nb=120`,
`max_iter_eq=40`, `tol_eq=2.5e-5`) on this macOS machine
(`.venv`: numpy 1.24.3 / numba 0.58.1 / miniconda MKL, AVX2 no AVX512).

## Reproduction band (Gate 0 outcome, revised to Gate 0')

The strict 1e-8 reproduction gate FAILED; the revised gate (contract
exactness + moments within 2% relative of the Torch record) PASSES:

- Contract exactness: theta, mechanism flags, and every override verified
  identical to the recorded M5 contract (see `gate0_diagnosis.json`).
- Local loss `8.998969403307` vs Torch `9.044422069071`: **~0.5% loss-level
  offset**.
- Per-moment offsets 1.1e-5 to 9.5e-3 in level, **up to ~1.7% relative**
  (largest: `prime30_55_childless_owner_minus_renter_mean_rooms`,
  `young_childless_renter_liquid_wealth_to_annual_gross_income_2535`).
- Both solves are strict (`residual <= tol_eq = 2.5e-5`): local `9.32e-06`
  vs Torch `1.48e-05`.
- The local solve is bit-identical across repeated runs on this machine
  (verified under two NUMBA thread settings), so the offset is a fixed
  cross-platform displacement, not noise.

## Cause hypothesis

Cross-platform floating-point non-reproducibility between this local
environment and the Torch environment that produced the record
(`module load anaconda3/2025.06`; different numpy/numba/MKL builds and SIMD
paths). The solver's `@njit(parallel=True)` reduction kernels and the loose
research tolerance of the equilibrium fixed point (`tol_eq=2.5e-5`, Brent
scalar refine) let build-level summation-order differences displace the
converged equilibrium slightly. Both points independently satisfy the same
convergence contract. The project's "bit-identical tight repeat" protocol has
only ever been validated same-machine; there is no precedent for 1e-8
cross-machine reproduction of a Torch record.

## Same-run-delta rule

Any local comparison across specifications (policy counterfactual vs.
benchmark, mechanism on vs. off) MUST compute both sides **in the same local
run** (same machine, same process, same environment), so the fixed platform
offset cancels in differences. Never mix a locally solved case with the
Torch-recorded baseline levels when computing deltas. Locally produced
LEVELS carry the ~0.5%/1.7% band above and must not be quoted as the
calibration's reported numbers; the Torch record remains the source for
levels.
