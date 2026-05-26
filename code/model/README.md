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
.venv/bin/python -m dt_cp_model.cli accounting-scale --setup fast --max-iter-eq 120 --quiet
.venv/bin/python -m dt_cp_model.cli scaled-equilibrium --setup fast --baseline-max-iter-eq 120 --max-iter-eq 80 --quiet
.venv/bin/python -m dt_cp_model.cli scaled-equilibrium --setup fast --baseline-max-iter-eq 1 --max-iter-eq 80 --outside-value -35.98609192692343 --force-full --quiet
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

Run the regression guard after changes touching entry, fertility, housing
demand, or equilibrium iteration:

```bash
.venv/bin/python tools/check_population_closure.py
```

See `PLAN.md` for the full implementation and optimization plan.
