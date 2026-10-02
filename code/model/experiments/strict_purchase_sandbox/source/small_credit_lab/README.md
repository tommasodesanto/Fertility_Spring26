# Refactor lab: stationary housing–fertility model (reviewed test package)

A readable, self-contained package for one fixed reference: **2007 stationary
reference — block0506, September 28 verified export** (checkpoint
repeat_0212, SHA `b15ba92d…`). It is a reviewed *test package*, not a new
baseline, not a calibration, and not an active-production replacement. The
dated transition and full reporting pipeline remain in the legacy code
(`code/model/tools`). Development history is archived in
`output/model/publication_refactor_20260929/`.

## What it solves

- **Model:** the maintained Markov-income sequential household model.
  One housing market, unit household mass, split birth-vintage adult entry,
  normalized population closure.
- **Stationary equilibrium:** with `tau_pay` fixed, the pension follows the
  analytic stationary rule. The single house price solves demand = static-elastic
  supply by direct Brent (`tol_eq = 2.5e-5`). The pension certificate is then
  checked on the solved distribution.
- **Inputs:** every primitive comes from the checkpoint through an explicit
  input bundle pinned by its SHA. There are no constructor defaults and no
  calibration; `psi_child = 0.1355551166583114` is fixed.
- **Credit:**
  - `--credit reference`: the checkpoint rule (renter age-taper rollover, DUE
    stayer and death floors).
  - `--credit corrected --d-bar D`: the reviewed upstream contract. Renters
    must save `b' ≥ −D`, with a zero floor at terminal or death-risk ages; an
    owner may sell into renting only if raw wealth after sale
    `b + (1−ψ)pH ≥ 0`; buyer and incumbent-owner rules are unchanged; borrowed
    negative principal is never rolled over.
  - **Warning:** with `D = 0`, two retained age-18 entrant cells
    (`b = −0.2558139535`) have no feasible choice. The run stops with the
    engine's census (exit 3); nothing is repaired, and corrected GE is refused.
    This needs an author decision on entrant debt.

## Running it (repo root, one core)

Python 3.13 **and** NumPy 2 are required to read the saved checkpoint pickle.
`requirements.txt` is a byte copy of
[`requirements_py313.freeze.txt`](../../../output/model/publication_refactor_20260929/local_env_v1/requirements_py313.freeze.txt); see
[`environment_py313_receipt.json`](../../../output/model/publication_refactor_20260929/local_env_v1/environment_py313_receipt.json) for the originating environment. Individual one-core runs locally, including overnight, are
permitted; Torch is for batches and parallel work.

```sh
export PY="$PWD/output/model/publication_refactor_20260929/local_env_v1/venv313/bin/python"
export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export PYTHONPATH="$PWD/code/model" MPLBACKEND=Agg
B="$PWD/output/model/publication_refactor_20260929/local_export_v1/inputs"
S=427e67a3d9dd663cd23c3f8533c55a1a64b4f9350d396c97b5c5bd4700bc90b7
OUT="$PWD/output/model/publication_refactor_20260929/local_runs/manual_ge_$(date +%Y%m%d_%H%M%S)"
"$PY" -m refactor_lab.run equilibrium --credit reference --initial-price-factor 1.05 \
  --reference-root "$PWD" --bundle "$B" --bundle-sha "$S" --out "$OUT"
```

Output (`receipt.json`, `ge/solution_arrays.npz`) is labelled
`preliminary_solve_only_not_certified`. The bundle is created once with
`python -m refactor_lab.export_inputs` from the pinned checkpoint.

## Model modules

| Path | Content |
|---|---|
| `engine/shared.py` | Primitives, income/entry helpers, closures |
| `engine/household.py` | Bellman stages, saving/tenure/fertility choices |
| `engine/distribution.py` | KFE, entry, statistics |
| `engine/equilibrium.py` | At-price solve and one-market GE |
| `engine/solver.py` | Facade re-exporting the four stages |
| `engine/kernels.py`, `parameters.py`, `utils.py`, `joint_nested.py`, … | Compiled kernels, credit floors, fertility nesting |
| `engine/e5f_stationary_paygo.py`, `e5f_social_security.py` | Pension rule and certificate |
| `inputs.py`, `export_inputs.py`, `credit.py`, `run.py` | Pinned inputs, credit binding, entry point |
| `verification/` | Oracles, benchmarks, provenance tools and drivers. `run.py` imports `verification.callcount`; instrumentation activates only with `--count-calls`. |

**Provenance.** The engine was extracted from the frozen source and the
separately reviewed three-file credit overlay; the materialized solver was
split without changes to its mathematical bodies. Reference mode preserves the
frozen borrowing behavior, while corrected mode applies the approved rule.
Indexed saving is the sole performance transformation after overlay
materialization. `exhaustive_saving_scalar` uses the indexed form with the same
candidates, order and strict tie rule. Its input was the
preserved scalar source (`19dceb70…`);
`verification/make_indexed_stage.py` produced `fe7d43af…`. The sequence and
promotion are recorded in `engine/transform_receipt.json`,
`engine/indexed_saving.diff` and `engine/promotion_receipt.json`. The promotion
is accepted for this fixed reference only.

## Verification (separate from the model)

`verification/` holds:
- the frozen-observer oracle (`acceptance_oracle.py`);
- same-machine baseline identity (`baseline_identity.py`);
- array comparison (`compare.py`);
- the budget supervisor (`budget_run.py`) and call counter (`callcount.py`);
- the native GE benchmark (`native_ge_benchmark.py`);
- provenance generators (`materialize.py`, `apply_split.py`, `make_indexed_stage.py`, `saving.py`);
- the phase driver (`verify.sh`, with the Torch wrapper `verify_torch.sh`);
- the source pins.

The certificate uses the frozen September 28 observer stack
(`run_fixed_price.py`, dated evaluation, gates, 14/31 tables, 17 standard
plots). That stack needs the intact reference source tree: locally, a
changed `tools/e5f_exact_policy_cache.py` makes it stop, deliberately and
without bypass. It therefore runs on Torch.

From the repository root, run the local suite or submit the certificate to
Torch. Set `LAB_SRC` to the staged source directory containing `refactor_lab/`;
the wrapper requires it because Slurm may relocate the submitted script.

```sh
export PY="$PWD/output/model/publication_refactor_20260929/local_env_v1/venv313/bin/python"
export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export PYTHONPATH="$PWD/code/model" MPLBACKEND=Agg
S=427e67a3d9dd663cd23c3f8533c55a1a64b4f9350d396c97b5c5bd4700bc90b7
V="$PWD/code/model/refactor_lab/verification/verify.sh"
VT="$PWD/code/model/refactor_lab/verification/verify_torch.sh"
MODE=local PHASE=acceptance BUNDLE_SHA="$S" bash "$V"

# Replace with the staged source directory containing refactor_lab/.
export LAB_SRC=/path/to/staged/source
sbatch --export=ALL,LAB_SRC,PHASE=fixed-price,BUNDLE_SHA="$S" "$VT"
```

## Evidence

- **Torch 18850069, fixed price.** Scalar and indexed engines, 2
  repetitions each: all 113 nested paths exact and finite, all 14 fit and 31
  parameter rows, 17 plot hashes exact.
- **Local native GE pair** (Apple M5 Pro, 1 thread, cold caches, identical
  instrumentation; same inputs, same start 1.05 × p_ref, both 4 household and
  5 KFE calls):
  - original engine 152.65 s, promoted engine 99.55 s of solve stage;
  - 90 arrays exactly equal, effective parameters equal, market and fiscal
    gates pass.
  - Excludes the historical calendar/table/plot certificate. It is a timing of
    identical work, not an economic result.
- **Renewal (diagnostic only).** GE relative gap 1.70e-6, versus 7.92e-7 in
  the reference `adult_entry_gate`. The 1e-6 threshold applies to the
  *difference* from the reference; it is not an absolute renewal gate, and
  the reference encodes none.
- **Torch full GE/reporting pair 18851943:** passed; 87 solution/shared
  arrays exact, 14/31 CSV rows byte-identical, and all 17 actual plot files
  hash-identical between engines. Market, fiscal and existing household/
  distribution gates pass. Solve stages: 275.95 s original, 204.87 s refactored.
  [Final report](../../../output/model/publication_refactor_20260929/REPORT.md)
  links the complete receipts, tables and plots.
- **Component suite (local):** 27 tests pass, including the provenance chain,
  credit contract and indexed-vs-original kernel on model columns.
