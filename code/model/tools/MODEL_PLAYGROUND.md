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

The runner uses one parameter file to perform one stationary GE solve. It searches for a price that clears the birth-renewal root; it does not search over parameter values. It then saves the complete case under `output/model/local_solution/cases/` and advances `latest` only after all storage and diagnostic outputs reopen successfully. Target scoring follows the GE solve. Each complete case contains 14 target-fit rows, 31 parameter rows, 17 native diagnostic figures, eight policy figures, and seven aggregate figures. The plotters load the validated selected case without solving. The public `model_run_io.load_run()` alias also resolves the production cache when called without a path; supply an explicit legacy run directory to read an older serialized format. The browser explorer requires the local server and the `explorer_cases.json` file from that case; see the production README for its command.

The separate [`best_params.py`](../parameters/best_params.py) and [`toy_params.py`](../parameters/toy_params.py) files hold complete inputs. Edit a copy such as `toy_params.py`: change `PARAMETERS["beta_annual"]` for discounting, or set `EXTERNAL_INPUTS["phi"] = [0.75] * 4` for a uniform financed share under the retained soft purchase rule. Its complementary 25% is the nominal non-financed share, not a strict liquid-cash threshold. The four entries must match; their period interpretation is not established. `unsecured_credit_limit` is renter debt capacity. Edit source earnings and payroll inputs rather than derived income or pension fields. Structural grid and entry-distribution edits are unsupported. Both files use the same production code and grid; they separate input values and saved outputs, not engine code.

From the project root, select the copy for the solve and for every later plot or explorer view:

```sh
"$PYTHON" "$PROJECT/code/model/run_model.py" --params toy_params.py
"$PYTHON" "$PROJECT/code/model/plot_model_policies.py" --params toy_params.py
"$PYTHON" "$PROJECT/code/model/plot_model_aggregates.py" --params toy_params.py
"$PROJECT/code/model/tools/start_model_explorer.command" --params toy_params.py
```

The exact canonical `code/model/parameters/best_params.py` routes to `output/model/local_solution/`. Most other files route to `output/model/experiments/<file-stem>/`; noncanonical files named `best_params.py` receive a source-specific directory. Plotters and the explorer read that file's latest completed case and do not run a model solve. If the case does not exist, they stop with a message; they do not show another parameter file's or a historical case by default. The full schema, validation, and export rules are in [`../parameters/README.md`](../parameters/README.md).

To open the fixed-price console directly with the toy inputs, pass the same file name to the playground:

```sh
"$PROJECT/code/model/tools/start_model_playground.command" --params toy_params.py
```

This can start before a stationary GE case exists. The console then reports that there is no saved result, leaves `sol=None`, and still creates `model` from `toy_params.py`. Run `model.solve(price=0.72)` for one fixed-price lifecycle solution. That does not find the renewal-price root or clear housing markets. Changes through `model.params` or `model.P` last only in that console; edit the parameter file to keep changes after exit.

The superseded fixed-price script is preserved as history under [`../../../calibration_archive/model_frontend_20261003/`](../../../calibration_archive/model_frontend_20261003/). The interactive interface below is the canonical companion for fixed-price partial-equilibrium inspection; `run_model.py` remains the stationary-GE entry point.

## Interactive fixed-price partial-equilibrium inspection

Launch [`start_model_playground.command`](start_model_playground.command) to open the canonical [`ModelPlayground`](model_playground.py) in the authenticated Python 3.13 environment. Initialization reads the selected input file and a matching saved result, if available, but performs no solve. If no matching result exists, it leaves `sol=None` and still creates `model` from the selected input file. It does not load a different or historical case automatically.

```python
model.show_parameters()
model.params["beta_annual"] = 0.98
changed = model.solve()
changed.aggregates()
changed.plot_policy(variable="consumption", age=30, income=4)
fig, axes = changed.plot_aggregates()
```

`model.solve()` runs one lifecycle solve at the current fixed price. It does not find a renewal-price root, fit calibration targets, or clear housing markets. Set a price directly with `model.solve(price=0.72)` or pass ordinary external or native inputs with `external_inputs={...}` or `native_overrides={...}`.

Use `model.params` for the ten displayed calibration coordinates. Supported direct `model.P` edits to non-parameter primitives, such as `sigma`, `R_gross`, `delta`, `tau_H`, and `H0`, are carried into the next solve; the equivalent explicit input dictionaries are also available. Direct edits to parameter-mapped fields, derived fields, entry distributions, structural grid fields, or `model.b_grid` fail explicitly. Edit source primitives instead of derived `income`, `pension`, `q`, or `user_cost_rate`. `model.solver` exposes the production equilibrium module. Use `model.reset()` to restore the selected input preset after an experiment.

`model.saved_result()` loads the hash-checked historical soft solution for an explicit comparison; it is not the current `latest` case. For example:

```python
baseline = model.saved_result()
changed = model.solve(overrides={"beta_annual": 0.98})
comparison = model.compare(baseline, changed)
comparison["overall"]
comparison["by_age"]
```

Result objects keep their own input and solution snapshots. `changed.g` is the realized post-tenure distribution, `changed.g_stay_distribution` is its owner-stayer subset, and `changed.g_beginning_distribution` is post-fertility and pre-tenure mass. Conditional policy arrays remain distinct from population aggregates. The saved 120-node wealth grid has not been certified as converged; see the [asset-grid diagnosis](../../../output/model/fixed_reference_economics_20260928/asset_grid_diagnosis_v1/README.md).
