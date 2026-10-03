# Local stationary production package

`code/model/production/` is the canonical package for the October 3, 2026 local stationary general-equilibrium workflow. The editable entry point is [`../run_model.py`](../run_model.py); it supplies explicit parameter and primitive-input dictionaries, then calls the shared workflow in `workflow.py`. The package contains the selected input snapshot and the stationary household, distribution, price, and equilibrium implementation, including `native_price.py` and `native_phase_b.py`. Frozen observer authentication and reporting remain read-only compatibility dependencies under `output/model/fixed_reference_economics_20260928/`; calibration also reads the authenticated chain-13 packet there. `engine_inventory.json` records copied engine files and source hashes.

The selected default is the author-adopted October 3 post-interest transaction-timing chain 13 input snapshot, with 120 wealth nodes and 9 income states. It retains the selected soft financing contract. The experimental `lowerA(m)` term is not enabled; hard and quarter purchase-rule variants are comparison cases, not adopted defaults. The old fixed-price interface is retained only as a historical snapshot under [`../../../calibration_archive/model_frontend_20261003/`](../../../calibration_archive/model_frontend_20261003/); `tools/model_run_io.py` is a compatibility loader: without a path it resolves the canonical production `latest` cache, while explicitly supplied legacy run paths use their saved format. The plotters read production storage directly.

## Run and inspect

Use the authenticated Python 3.13 environment, which is compatible with the retained reference pickle inputs:

```sh
PROJECT=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
PYTHON="$PROJECT/output/model/publication_refactor_20260929/local_env_v1/venv313/bin/python"
"$PYTHON" "$PROJECT/code/model/run_model.py"
```

The runner calls one stationary GE workflow. Its price root is the birth-renewal residual
\(B_{adj}/(2.1 E)-1\); target scoring occurs after the solve and is not part of the root. With the default `fixed_h0` closure, the physical supply coefficient remains fixed. Let \(D(q)\) be normalized housing demand, \(u_c\) the user cost, and \(\xi\) the supply elasticity. The implied population-one supply coefficient is \(H_0^{\mathrm{implied}}=D(q)/(u_c q/\bar r)^\xi\), and the fixed-supply population scale is \(N=H_0^{\mathrm{fixed}}/H_0^{\mathrm{implied}}\). Both are calculated from the same selected policies and price, without a second solve; this accounting relies on conditional scale independence of policies and the fixed-payroll closure. The `population_one` closure is the production calibration normalization and is also used for comparison; it fixes normalized household population at one and derives the housing coefficient at the solved price.

Each attempt receives a unique folder under `output/model/local_solution/cases/`. A successful case contains the native solver report, a complete target-fit table, a parameter table, cached native arrays and metadata, 17 native diagnostic figures, 8 policy figures, 7 aggregate figures, and browser-explorer assets. The `latest` pointer advances only after the result and artifacts reopen successfully. Failed attempts remain in their case folder and do not replace `latest`.

The existing plotters load `output/model/local_solution/latest` without solving again:

```sh
"$PYTHON" "$PROJECT/code/model/plot_model_policies.py"
"$PYTHON" "$PROJECT/code/model/plot_model_aggregates.py"
```

The browser explorer needs its local server; opening its HTML directly is insufficient. From the project root, run `"$PYTHON" code/model/tools/economics_explorer.py --config output/model/local_solution/latest/explorer_cases.json` and use the local address it prints.

Edit the explicit `PARAMETERS`, `EXTERNAL_INPUTS`, and `NATIVE_OVERRIDES` dictionaries in `run_model.py` for a local experiment. The gross-earnings controls are `w_hat`, `income_age_profile`, and `tau_pay`; the accepted fixed-payroll mapping derives disposable income and pension income from them. Do not edit derived income or pension fields directly. Period-native `R_gross`, `delta`, and `tau_H` are authoritative; do not supply duplicate annual versions. Structural edits to entry distributions or grid dimensions are unsupported by this saved input snapshot and fail explicitly. Other valid edits create a new run and do not alter the target contract or paper baseline.

## Calibration entry

`calibration.py` uses the same stationary GE function, validates the complete pinned target/weight fingerprint before solving, and scores the 10 informative moments after each GE solve. Its search coordinates and acceptance bounds are derived from the authenticated chain-13 packet. The 14-row full target-fit table and 31-row parameter table are written for inspection. Use the owning calibration driver and its launch contract; this package README does not authorize a new search.

## Verification scope and limits

The local deployment passed its October 3 verification within the documented scope. Against a fresh execution of the authenticated original reference, both the unchanged chain-13 case and a case lowering annual \(\beta\) by 0.001 match all 91 arrays, all fields of the 14-row target-fit table, all 31 numeric parameter fields and bounds, and all 17 standard-plot hashes. The comparison receipts list 13 enumerated descriptive-only role/status differences and no numeric differences: [unchanged](../../../output/model/production_deployment_20261003/compare_unchanged/comparison.json) and [beta change](../../../output/model/production_deployment_20261003/compare_beta/comparison.json). Six synthetic workflow checks and a 245-field default/repeated-context check passed without lifecycle solves. The beta case is a numerical parity check, not an estimate or adoption.

The default `run_model.py` case is saved at `output/model/local_solution/latest`; it contains the 17 standard figures, 8 policy figures, 7 aggregate figures, and explorer configuration. The cached plotters and the live explorer's four HTTP routes were checked with zero model solves. An external \(\sigma=2.01\) GE input experiment propagated to 58 changed arrays and a changed price; it is not recalibration or adoption ([receipt](../../../output/model/production_deployment_20261003/external_input_propagation.json)). The algebraic scale exercise reused the beta case's policies and price without a second solve: \(H_0^{\mathrm{implied}}=6.39859337\), versus fixed \(H_0=6.405693596\), giving \(N=1.00110965\). Its normalized distributions still sum to one; absolute totals scale by \(N\) ([receipt](../../../output/model/production_deployment_20261003/scale_exercises.json)). This relies on conditional scale independence of policies and fixed-payroll closure.

Four maintained transition modules import successfully, and the 336-file live cluster bundle and independent legacy driver remain hash-checked and unchanged. No dynamic transition was replayed. Superseded frontends alone were archived; authenticated observer and oracle bundles remain read-only production dependencies. This verification establishes the local deployment within the checks above. It does not establish a new calibration, alter the chain-13 working anchor, or certify a dated transition, global optimum, grid adequacy, or paper baseline. See [the live calibration status](../../../CALIBRATION_STATUS.md) for unresolved research questions.
