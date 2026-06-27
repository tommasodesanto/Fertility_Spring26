# Python Model Codebase

This folder contains the active Python implementation of the current
discrete-time center-periphery fertility model. The former MATLAB code is now
archived for historical reference and parity checks.

Use the local virtualenv created for the port:

```bash
cd code/model
.venv/bin/python -m dt_cp_model.cli smoke
```

The virtualenv uses conda's NumPy and has `numba==0.58.1` installed for the
compiled golden-section savings kernels. If the environment has to be rebuilt:

```bash
/Users/tommasodesanto/miniconda3/bin/python -m venv --system-site-packages .venv
.venv/bin/python -m pip install 'numba==0.58.1'
```

Useful commands:

```bash
.venv/bin/python -m dt_cp_model.cli smoke --quiet
.venv/bin/python -m dt_cp_model.cli solve-theta --setup fast --theta x0 --max-iter-eq 120 --quiet
.venv/bin/python -m dt_cp_model.cli benchmark --setup fast --max-iter-eq 1 --force-full --quiet
.venv/bin/python -m dt_cp_model.cli objective --setup fast --max-iter-eq 1 --inv-iters 1 --quiet
.venv/bin/python tools/compare_benchmarks.py benchmarks/fast_ge_iter1.json benchmarks/matlab_fast_ge_iter1.json
.venv/bin/python tools/compare_benchmarks.py benchmarks/solve_theta_x0_fast_compiled_forward.json benchmarks/matlab_solve_theta_x0_fast.json
.venv/bin/python tools/check_population_closure.py
.venv/bin/python tools/export_room_pattern_diagnostics.py --case hR8=../../output/model/reduced_target_overnight_20260527/records/hR8_default_best.json,8.0 --case hR6=../../output/model/reduced_target_overnight_20260527/records/hR6_micro_best.json,6.0
.venv/bin/python tools/run_policy_counterfactuals_from_record.py --outdir ../../output/model/policy_counterfactuals_live_20260528_hR8_default
.venv/bin/python -m dt_cp_model.cli accounting-scale --setup fast --max-iter-eq 120 --quiet
.venv/bin/python -m dt_cp_model.cli scaled-equilibrium --setup fast --baseline-max-iter-eq 120 --max-iter-eq 80 --quiet
.venv/bin/python -m dt_cp_model.cli scaled-equilibrium --setup fast --baseline-max-iter-eq 1 --max-iter-eq 80 --outside-value -35.98609192692343 --force-full --quiet
```

## Fast Intergen One-Run Review

For the one-market intergenerational strand under
`intergen_housing_fertility/`, the default interactive inspection path should
stay small: one candidate, one solve or trusted solution cache, one quick plot
packet. Do not add `--run-policy-cases` unless the goal is explicitly a GE
counterfactual exercise.

From this directory:

```bash
PYTHONPATH=$PWD NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python tools/build_intergen_mechanics_packet.py \
  --source ../../output/model/intergen_current_review/source_candidate.json \
  --outdir ../../output/model/intergen_current_review/quick \
  --target-set candidate_replacement_roomgap_14moment_tfr192_v1 \
  --J 16 --Nb 60 --income-states 5 --n-house 5 --hR-max 6.0 \
  --max-iter-eq 10 --interp-method linear \
  --clean-outdir --no-csv \
  --skip-standard-diagnostics --skip-contact-sheet --quick-first-look-only
```

If a trusted `solution_cache.pkl` already exists, rebuild the quick plots
without re-solving:

```bash
PYTHONPATH=$PWD NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python tools/build_intergen_mechanics_packet.py \
  --source ../../output/model/intergen_current_review/source_candidate.json \
  --outdir ../../output/model/intergen_current_review/quick \
  --target-set candidate_replacement_roomgap_14moment_tfr192_v1 \
  --J 16 --Nb 60 --income-states 5 --n-house 5 --hR-max 6.0 \
  --max-iter-eq 10 --interp-method linear \
  --solution-cache ../../output/model/intergen_current_review/solution_cache.pkl \
  --refresh-plots-from-cache \
  --clean-outdir --no-csv \
  --skip-standard-diagnostics --skip-contact-sheet --quick-first-look-only
```

With `--no-csv`, the quick packet writes only the target table, moments,
solution summary, liquid- and total-wealth first-look policy figures, and
liquid- and total-wealth density figures. Use the full packet only for a slower
audit.

To audit the entry-wealth atom and owner-rung total-wealth jaggedness from the
current point-entry and five-node-entry caches:

```bash
PYTHONPATH=$PWD .venv/bin/python tools/audit_intergen_entry_atom.py \
  --outdir ../../output/model/intergen_current_review/atom_audit
```

The intended model boundary is:

```python
from dt_cp_model import solve_theta

sol, P, p_eq = solve_theta(theta, setup_mode="fast")
```

The SMM objective should remain a wrapper around this same parameter-vector
application path, not a separate interpretation of `theta`.

## Population Closure

The live calibration uses the benchmark-normalized outside-option closure:
`P.population_closure = "outside_option_benchmark_normalized"`.

In the benchmark, the model solves for a normalized stationary distribution,
calibrates the outside value to a target city-entry probability \(q^E\), and
sets the outside-born potential entrant mass mechanically so that benchmark
scale is \(S=1\):
\[
S E_0(p)=q^E(p)\left[M+S B_0(p)\right],
\qquad
M=\frac{E_0(p)}{q^E(p)}-B_0(p).
\]

Counterfactuals should hold the benchmark \(M\), \(\bar W^E\), and entrant
taste scale fixed, then let the implied stationary scale move with prices and
fertility. The older `normalized`, `renewal_valve`, and
`accounting_scale_prices` closures remain in code for diagnostics and
comparison, but they are not the current benchmark calibration closure.

The policy counterfactual runner follows that boundary: it first re-solves the
record under `outside_option_benchmark_normalized`, recovers the benchmark
outside objects, and then solves policy cases with
`population_closure="accounting_scale_prices"` while holding those outside
objects fixed. By default it runs the two current housing-supply cases:
all-location `H0 +10%` and center `H0 +10%`. Add repeated `--case ...`
arguments to run credit, tax, or transfer cases from the same driver.

Run the regression guard after changes touching entry, fertility, housing
demand, or equilibrium iteration:

```bash
.venv/bin/python tools/check_population_closure.py
```

See `PLAN.md` for the full implementation and optimization plan.
