# Surrogate-assisted calibration prototype (one-market intergen strand)

Status: **diagnostic methods prototype**, not a production SMM calibration. It
does not change the model, the target system, or any calibration state. It is an
alternative *outer-loop* for the existing one-market intergenerational
housing-fertility model in `code/model/intergen_housing_fertility/`.

## Why

The model solve is cheap (~40s) and low-dimensional; the expensive, ill-behaved
part is the SMM outer loop — blind global differential evolution on a cluster
plus a finite-difference Jacobian over a near-discrete objective (the June-18
identification ledger reports condition number ~2.7e4 at the core point and
local rank deficiency at the room-cost point). This prototype attacks that outer
loop with a **multi-output moment emulator**:

- a Gaussian process per target moment, `theta -> m(theta)`, trained on the
  existing 15k evaluation records;
- from which we read a **smooth** Jacobian `d m / d theta` analytically (the
  ledger's "Immediate Next Check"), comparable head-to-head with the FD audit;
- and run **Bayesian optimization** (expected improvement) to propose the next
  batch of model evaluations, plus a Pareto-frontier probe of whether the
  ledger's conflicting moments are jointly reachable.

Two honest limits, stated up front: a smooth surrogate can paper over genuine
non-identification, so the surrogate Jacobian *complements* the structural audit
rather than replacing it; and no outer loop repairs an incoherent target system
— if the frontier moments are jointly unreachable, faster search just converges
to the same wall. The frontier probe is designed to expose exactly that.

## What's here

| File | Role |
|---|---|
| `gp.py` | Pure-numpy ARD (anisotropic-RBF) GP. Analytic LML gradients, Adam fit, analytic posterior-mean gradient. Has a self-test (`python -m intergen_surrogate_calibration.gp`). |
| `data.py` | Loads/curates the 15k Tier-1 records. Imports bounds, target sets, weights, and `diagnostic_loss` from the model package (no re-encoding). |
| `emulator.py` | `MomentEmulator`: one GP per moment + emulated loss + MC expected-loss + K-fold cross-validation. |
| `identification.py` | Surrogate Jacobian with the FD audit's exact target-normalization; SVD rank/condition/collinearity; ARD relevance; head-to-head with `intergen_sensitivity_jacobian_20260618`. |
| `bo.py` | MC-EI acquisition, diverse batch proposals, exploit point, and a Pareto frontier over conflicting moments. |
| `run_surrogate_calibration.py` | Driver: trains the emulator, runs CV + identification + BO, writes a packet + plots + `REPORT.md`. |
| `validate_proposals.py` | Solves proposals and reproduction-checks records with the **real** `run_model_cp_dt`. |

## No new dependencies

The model venv ships numpy/pandas/matplotlib/numba but **not** scipy or
scikit-learn, so the GP is hand-rolled in numpy and everything runs in the model
`.venv`. The real-solver validation runs in the same env.

## Run

```bash
cd code/model

# 1) GP unit self-test (gradients vs finite differences, OOS recovery)
.venv/bin/python -m intergen_surrogate_calibration.gp

# 2) Full packet: emulator + cross-validation + identification + BO
.venv/bin/python -m intergen_surrogate_calibration.run_surrogate_calibration \
    --target-set candidate_no_timing_v0 --n-keep 1500 --restarts 2 --iter 80

# 3) Validate proposals against the real solver (~40-50s per solve)
.venv/bin/python -m intergen_surrogate_calibration.validate_proposals \
    --target-set candidate_no_timing_v0 --limit 4 --repro 2
```

Outputs: `output/model/intergen_surrogate_calibration/<target_set>/`
(`REPORT.md`, `VALIDATION.md`, JSON records, and PNG diagnostics).

## Data sources (Tier-1, full theta recoverable)

- `output/model/cluster_pulls/results_intergen_housing_fertility_intergen_fast_globalde_hown5_refinement_v1_2h_20260617/` (7,901)
- `output/model/cluster_pulls/intergen_replacement_cluster_wave_20260618/` (7,200)

Records store the raw moment dict, so the emulator predicts moments and any
registered target set's loss is recomputed with the model's own
`diagnostic_loss`. `candidate_no_timing_v0` is the primary demonstration because
all 13 of its moments are logged in ~15k records and it is the exact target
system used by the finite-difference Jacobian audit.

## Path to scale-up (not done here)

- Inducing-point / sparse GP to use all 15k points instead of a curated ~1.5k.
- Active-learning loop: propose batch → solve on Torch → append → refit.
- Differentiable solver (autodiff) for an *exact* Jacobian if the emulator
  Jacobian and FD audit ever disagree materially.
