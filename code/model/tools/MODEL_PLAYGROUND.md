# Run and inspect the current production model

The current editable runner is [`../run_model.py`](../run_model.py). It calls the canonical stationary-GE workflow in [`../production/README.md`](../production/README.md) using the author-adopted October 3 post-interest chain-13 input snapshot. Use the Python 3.13 interpreter shown below; the older project venv may not read the reference pickle inputs.

```sh
PROJECT=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
PYTHON="$PROJECT/output/model/publication_refactor_20260929/local_env_v1/venv313/bin/python"
"$PYTHON" "$PROJECT/code/model/run_model.py"
"$PYTHON" "$PROJECT/code/model/plot_model_policies.py"
"$PYTHON" "$PROJECT/code/model/plot_model_aggregates.py"
```

The runner performs one birth-renewal GE root, then saves the complete case under `output/model/local_solution/cases/` and advances `latest` only after all storage and diagnostic outputs reopen successfully. Target scoring follows the GE solve. Each complete case contains 14 target-fit rows, 31 parameter rows, 17 native diagnostic figures, eight policy figures, and seven aggregate figures. The plotters load the validated latest case without solving. The public `model_run_io.load_run()` alias also resolves the production cache when called without a path; supply an explicit legacy run directory to read an older serialized format. The browser explorer requires the local server and the `explorer_cases.json` file from that case; see the production README for its command.

Edit `PARAMETERS`, `EXTERNAL_INPUTS`, and `NATIVE_OVERRIDES` at the top of `run_model.py` to make a local experiment. The gross-earnings controls are `w_hat`, `income_age_profile`, and `tau_pay`; the accepted fixed-payroll mapping derives disposable income and pension income. Do not edit derived income or pension fields directly. Period-native `R_gross`, `delta`, and `tau_H` are authoritative; do not enter duplicate annual values. Structural entry-distribution or grid-dimension edits are unsupported and fail explicitly. Valid edits create a new run; they do not alter calibration targets or paper baselines. The local workflow passed its October 3 verification within the tested scope. Exact same-host parity against the executed reference covers the unchanged case and a \(\beta\)-minus-0.001 case: all 91 arrays, the complete 14-row target table, numeric values and bounds in all 31 parameter rows, and 17 standard-plot hashes match; only enumerated descriptive metadata differs. The verified `latest` cache contains 17 standard, 8 policy, and 7 aggregate figures. Cached plots and four live explorer routes were checked without solves. See [the comparison receipts](../../../output/model/production_deployment_20261003/compare_unchanged/comparison.json) and [production guide](../production/README.md). This does not establish a new calibration or dynamic transition.

The superseded fixed-price script is preserved as history under [`../../../calibration_archive/model_frontend_20261003/`](../../../calibration_archive/model_frontend_20261003/). The interactive interface below is the canonical companion for fixed-price partial-equilibrium inspection; `run_model.py` remains the stationary-GE entry point.

## Interactive fixed-price partial-equilibrium inspection

Launch [`start_model_playground.command`](start_model_playground.command) to open the canonical [`ModelPlayground`](model_playground.py) in the authenticated Python 3.13 environment. Initialization reads production inputs and the saved result but performs no solve. The global `sol` uses `output/model/local_solution/latest` when available; otherwise it loads the hash-checked historical soft case. The public `ModelPlayground` class always uses the current production inputs and engine.

```python
model.show_parameters()
model.params["beta_annual"] = 0.98
changed = model.solve()
changed.aggregates()
changed.plot_policy(variable="consumption", age=30, income=4)
fig, axes = changed.plot_aggregates()
```

`model.solve()` runs one lifecycle solve at the current fixed price. It does not find a renewal-price root, fit calibration targets, or clear housing markets. Set a price directly with `model.solve(price=0.72)` or pass ordinary external or native inputs with `external_inputs={...}` or `native_overrides={...}`.

Use `model.params` for the ten displayed calibration coordinates. Supported direct `model.P` edits to non-parameter primitives, such as `sigma`, `R_gross`, `delta`, `tau_H`, and `H0`, are carried into the next solve; the equivalent explicit input dictionaries are also available. Direct edits to parameter-mapped fields, derived fields, entry distributions, structural grid fields, or `model.b_grid` fail explicitly. Edit source primitives instead of derived `income`, `pension`, `q`, or `user_cost_rate`. `model.solver` exposes the production equilibrium module. Use `model.reset()` to restore defaults after an experiment.

`model.saved_result()` loads the hash-checked historical soft solution for an explicit comparison; it is not the current `latest` case. For example:

```python
baseline = model.saved_result()
changed = model.solve(overrides={"beta_annual": 0.98})
comparison = model.compare(baseline, changed)
comparison["overall"]
comparison["by_age"]
```

Result objects keep their own input and solution snapshots. `changed.g` is the realized post-tenure distribution, `changed.g_stay_distribution` is its owner-stayer subset, and `changed.g_beginning_distribution` is post-fertility and pre-tenure mass. Conditional policy arrays remain distinct from population aggregates. The saved 120-node wealth grid has not been certified as converged; see the [asset-grid diagnosis](../../../output/model/fixed_reference_economics_20260928/asset_grid_diagnosis_v1/README.md).
