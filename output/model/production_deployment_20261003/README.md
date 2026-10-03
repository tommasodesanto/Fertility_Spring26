# October 3 stationary production deployment: verified readout

## Verified facts and source identity

- **Reference:** adopted post-interest transaction timing, soft financing, chain 13. The exact ten coordinates and native price start `0.7760569760205563` come from the [selected run receipt](../fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/collection/production_alternative_chain_13/run/completed.json) and its [passed native postcheck](../fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/collection/production_alternative_chain_13/run/native_postcheck/completed.json).
- **Executed primitives:** 120 wealth nodes, 9 income states, four-year periods, housing measured in rooms; provisional nonnegative-mean entrant wealth mapping with mean wealth `0.18651967924681834`. [Local input snapshot and source hashes](../../../code/model/production/reference_inputs/bundle.json).
- **Matched numerical checks:** unchanged and lower-beta production results each match their fresh authenticated reference exactly in 91 saved solution/shared arrays, all 14 target rows, all numeric fields/restrictions in 31 parameter rows, closure residuals, and all 17 standard diagnostic PNG hashes. [Unchanged comparison](compare_unchanged/comparison.json), [beta comparison](compare_beta/comparison.json).
- **Runtime:** separate fresh processes on the same macOS ARM64 host; Python 3.13.15, NumPy 2.2.6, SciPy 1.15.3, Numba 0.61.2 and Matplotlib 3.10.3. Numba, OpenMP, OpenBLAS and MKL were each limited to one thread. Full runtime identity is recorded in each execution receipt.
- **Canonical sources:** [production package](../../../code/model/production/README.md), [18-file engine hash inventory](../../../code/model/production/engine_inventory.json), [native price/root source inventory](../../../code/model/production/native_core_inventory.json), and [source/retention map](../../../code/model/production/relocation_manifest.json). The selected snapshot pins the complete target fingerprint `db60605ef444b747b2ed7f9482b6b8c90ef5e8b1d75c49a3ac1ff6f4275c9ba1` and weight fingerprint `2391cd2d4a39a6669a405be34ff14116ad314354353173cb619f7fa7c66043b0`.

This verifies the deployed stationary workflow at the tested inputs. It does not certify a global calibration optimum, grid adequacy, a dated transition, or publication readiness. The working scientific qualifications remain in [CALIBRATION_STATUS.md](../../../CALIBRATION_STATUS.md).

## Fresh stationary executions and exact comparisons

The population-one normalization fixes household scale at \(N=1\), roots the birth-renewal residual, and derives the housing supply coefficient. Both backends used this normalization in the matched comparisons. The changed-beta case lowers annual beta from `0.9663191380998087` to `0.9653191380998087`; it changes no other supplied primitive.

| Case | Reference / production lifecycle solves | Price | Birth-renewal residual | Housing residual | PAYGO residual | Result |
|---|---:|---:|---:|---:|---:|---|
| Unchanged | 2 / 2 | 0.7760569760205563 | −6.82106593430376e−11 | 0 | 2.846866047437337e−14 | [Exact numerical parity](compare_unchanged/comparison.json) |
| Annual beta lower by 0.001 | 7 / 7 | 0.7759400687622401 | −3.154748429157195e−8 | 0 | 2.9518734016460887e−14 | [Exact numerical parity](compare_beta/comparison.json) |
| Risk aversion 2.01 | Production: 9 | 0.778424953568571 | −9.579337323373238e−11 | 0 | 2.870201015039283e−14 | [External-input propagation](external_input_propagation.json) |

Full economic readouts are linked below; no target or estimated-parameter row is omitted from those tables.

| Case | Complete 14-row target fit | Complete 31-row parameter table | Execution receipt |
|---|---|---|---|
| Unchanged reference | [Target, model, gap, weight, loss](reference_unchanged/report/target_fit.csv) | [Estimates and restrictions](reference_unchanged/report/parameters.csv) | [Receipt](reference_unchanged/receipt.json) |
| Unchanged production | [Target, model, gap, weight, loss](production_unchanged_v2/report/target_fit.csv) | [Estimates and advisory reference bounds](production_unchanged_v2/report/parameters.csv) | [Receipt](production_unchanged_v2/receipt.json) |
| Lower-beta reference | [Full fit](reference_beta/report/target_fit.csv) | [Full parameters](reference_beta/report/parameters.csv) | [Receipt](reference_beta/receipt.json) |
| Lower-beta production | [Full fit](production_beta/report/target_fit.csv) | [Full parameters](production_beta/report/parameters.csv) | [Receipt](production_beta/receipt.json) |
| Risk-aversion experiment | [Full fit](production_external/report/target_fit.csv) | [Full parameters](production_external/report/parameters.csv) | [Receipt](production_external/receipt.json) |

Each matched comparison checks the canonical archive `selected_repeat/stage/solution_arrays.npz`: 91 arrays, including 20 shared arrays. It requires nonempty inventories and named value, consumption, saving, housing, tenure, fertility and distribution arrays. It also checks every table column, all common closure fields and actual hashes for the unchanged 17-plot set. Missing or unequal required artifacts fail the CLI.

The parameter CSVs are **not byte-identical**. Each comparison records 13 enumerated descriptive differences: ten supplied-primitive roles, the derived housing-supply role, the derived child-benefit coefficient role, and the corresponding closure description. All estimates, reference estimates, bounds, near-bound flags and fixed restrictions match. All target-table metadata and numeric values match; the target CSV bytes also match. Every permitted before/after description appears in the comparison receipt. Unknown differences remain fatal. Additional scale metadata, if present, must satisfy the explicitly checked algebraic identities.

The risk-aversion experiment confirms that `sigma=2.01` reaches the executed parameter report and changes 58 saved arrays, including policies and the distribution. The price changes from `0.7760569760205563` to `0.778424953568571`. This is an ordinary input experiment, not recalibration or adoption. There is no matched reference external-input execution: its authenticated expected-parameter observer would need a separately scoped extension.

## Editable runner, saved results and inspection

The final [run_model.py](../../../code/model/run_model.py) execution used the default fixed-\(H_0\) normalization and completed at the selected price. Its validated cache is [output/model/local_solution/latest](../local_solution/latest/SUMMARY.md), pointing to case `20261003T175652812716Z_b1c72f13`. Housing residual is zero, the birth-renewal residual is `−6.82106593430376e−11`, and the PAYGO residual is `2.846866047437337e−14`.

The successful cache contains 17 standard diagnostic figures, 8 policy figures and 7 aggregate figures, plus the [complete target fit](../local_solution/latest/target_fit.csv) and [parameter table](../local_solution/latest/parameters.csv). Both standalone plotter commands were rerun from the saved cache with zero model solves. The [HTTP explorer check](explorer_check.json) obtained status 200 for all four tested routes: `/`, `/api/meta`, a policy slice, and aggregate data.

The [cache metadata correction](cache_metadata_correction.json) fixed an inaccurate statement that the baseline fiscal mapping had changed. The supplied and effective baseline fiscal/entry primitives equal the authenticated reference. The numerical archive hash and `latest` pointer remained unchanged; this correction required no solve or replot.

One operator mistake is recorded for auditability: before argument parsing was added, I invoked `run_model.py --help`, and the runner ignored the flag and completed a second native solve. The extra case is retained at [20261003T181552418704Z_511ca917](../local_solution/cases/20261003T181552418704Z_511ca917/); the intended published case remains [20261003T175652812716Z_b1c72f13](../local_solution/cases/20261003T175652812716Z_b1c72f13/), and `latest` points to it. Both `native_result.npz` files have SHA-256 `27cae50d61fe07ce44733360b9f1583b4e229581661694a39e9b579820d5e8b6`, so the retained extra is byte-identical to the prior result. Its observer receipt records two lifecycle solves, about 27 seconds from case reservation through metadata creation, and the runner caps Numba, BLAS and OpenMP at one thread. The CLI now handles `--help` before importing the workflow and rejects unknown arguments; a pure test confirms neither path creates an output case.

The zero-solve [interactive input check](interactive_input_check.json) confirms direct `P.sigma` and `P.R_gross` edits reach the canonical adapter, structural grid edits fail explicitly, and package-relative plotting imports work. Seven storage/workflow tests passed, and the calibration adapter self-test passed for supplied-primitive/grid copying, the common solver call, serializable optimizer fields and pre-solve target fingerprint rejection. A pure input check using the actual runner dictionaries also retained a 1% gross-wage edit and the existing entrant wealth distribution, deriving disposable income and pensions under the explicitly accepted fixed-payroll mapping. This was an input-mapping check, not a numerical wage experiment.

## Ordinary inputs, calibration bounds and scale accounting

[Input validation](input_validation.json) establishes that annual beta `0.995` and parent housing floor `3.0` are preserved by the ordinary input path even though they exceed the reference calibration upper bounds `0.99` and `2.6`. It used zero model solves and does not establish equilibrium convergence at those values. Nonuniform financed shares are rejected because the active native contract does not support them.

The [input binder](../../../code/model/production/inputs.py), [stationary observer](../../../code/model/production/native_phase_b.py) and [calibration adapter](../../../code/model/production/calibration.py) keep the distinction explicit: reference bounds and near-bound labels are advisory reporting metadata in ordinary GE. Calibration enforces its search bounds and derived-\(H_0\) acceptance interval `[0.2,80]`; GE does not enforce those optimizer bounds. Both use the same stationary solver. The current alternative-timing calibration driver defaults to production, while `--historical-reference` selects the retained authenticated evaluator. Native postcheck subprocesses inherit the route.

The [saved-array scale exercise](scale_exercises.json) uses changed-beta stationary policies and the normalized distribution. With retained \(H_0=6.40569359569417\) and implied population-one coefficient `6.398593372825029`, algebraic rescaling gives \(N=1.0011096537090942\), zero housing residual and unchanged relative birth-renewal residual `−3.154748429157195e−8`. It is an algebraic check with zero model solves, not another GE solve or a fresh forward-distribution replay. Policies are reused unchanged and absolute household quantities scale with \(N\).

## Preserved dependencies, failed attempts and limits

The household solver, forward distribution and stationary price/GE code execute from `code/model/production/`. Frozen observer authentication, empirical measurements and accounting certificates remain read-only compatibility dependencies under historical output folders; the deployment does not eliminate them. Source identities and retained locations are recorded in the inventories above.

The [deployment-preservation receipt](live_deployment_preservation.json) verifies the immutable staged legacy deployment at commit `537eaae8`: original and continuation drivers still match its inventory, its 336-file staged archive contains the legacy evaluator, and canonical production was not inserted. Existing pinned bundles and comparison sources remain preserved. Four maintained transition entrypoints passed [import checks](transition_imports.json) with zero solves. That establishes import compatibility only; no transition replay or transition acceptance claim follows.

Earlier attempts remain inspectable: [production_unchanged_import_failed/](production_unchanged_import_failed/) failed before solving because direct-script execution lacked package-relative imports; [production_unchanged/](production_unchanged/) failed after one lifecycle solve at an observer object-identity check. The corrected source callback preserves the existing economic gates. The final successful unchanged run is `production_unchanged_v2`, followed by the successful beta and external-input runs. These failures were retained, not overwritten.

## Commands

Use the project root and the retained Python environment. Importing the runner or package is inert. Each verification execution starts one fresh process, with a 1,800-second budget and 32-lifecycle cap. The commands below identify the completed runs; existing output directories are refused, so choose fresh output names for a new replay.

```sh
export PYTHONPATH=code/model
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
MODEL_PY=output/model/publication_refactor_20260929/local_env_v1/venv313/bin/python

"$MODEL_PY" -m production.verification execute --backend reference --case unchanged --out output/model/production_deployment_20261003/reference_unchanged --budget-seconds 1800 --max-lifecycle 32
"$MODEL_PY" -m production.verification execute --backend production --case unchanged --out output/model/production_deployment_20261003/production_unchanged_v2 --budget-seconds 1800 --max-lifecycle 32
"$MODEL_PY" -m production.verification compare --left output/model/production_deployment_20261003/reference_unchanged --right output/model/production_deployment_20261003/production_unchanged_v2 --out output/model/production_deployment_20261003/compare_unchanged

"$MODEL_PY" -m production.verification execute --backend reference --case beta --out output/model/production_deployment_20261003/reference_beta --budget-seconds 1800 --max-lifecycle 32
"$MODEL_PY" -m production.verification execute --backend production --case beta --out output/model/production_deployment_20261003/production_beta --budget-seconds 1800 --max-lifecycle 32
"$MODEL_PY" -m production.verification compare --left output/model/production_deployment_20261003/reference_beta --right output/model/production_deployment_20261003/production_beta --out output/model/production_deployment_20261003/compare_beta

"$MODEL_PY" -m production.verification execute --backend production --case external --out output/model/production_deployment_20261003/production_external --budget-seconds 1800 --max-lifecycle 32
"$MODEL_PY" -m production.calibration --self-test
"$MODEL_PY" -m pytest -q code/model/production/test_storage.py

"$MODEL_PY" code/model/run_model.py
"$MODEL_PY" code/model/plot_model_policies.py
"$MODEL_PY" code/model/plot_model_aggregates.py
"$MODEL_PY" code/model/tools/economics_explorer.py --config output/model/local_solution/latest/explorer_cases.json
```

The plotter and explorer commands inspect saved results. `run_model.py` launches a new stationary solve. No new cluster job, calibration search, or dated transition was launched for this deployment.

## External parameter files and personal experiments

The runner now selects `code/model/parameters/best_params.py` or an independent
editable copy such as `toy_params.py`. The [current code map](../../../code/model/README.md)
and [parameter guide](../../../code/model/parameters/README.md) give the commands.
The default still uses the same verified chain-13 inputs. Other files have
separate latest pointers, plots and explorer configurations; edited presets are
saved with their exact text and hash, and cache readers disclose later edits.

The [workflow receipt](parameter_file_workflow.json) records 40 passing focused
tests, a zero-solve canonical calibration initialization, both cached plotter
commands and the toy browser identity. One complete experimental beta-minus-0.001
GE run converged with 17 standard, 8 policy and 7 aggregate plots. Its full
[target table](../experiments/toy_params/latest/target_fit.csv),
[parameter table](../experiments/toy_params/latest/parameters.csv) and
[summary](../experiments/toy_params/latest/SUMMARY.md) remain separate from the
production cache; the latter's original latest pointer and NPZ hash were unchanged.
This is a parameter experiment, not a new calibration or accepted specification.

The initial end-to-end attempts exposed frozen authentication dependencies on
archived source and a native verification reserve exceeding the proposed
600-second cap. All 24 cleanup-affected frozen pins now match again; the
[archive record](../../../calibration_archive/model_legacy_20261003/README.md)
identifies compatibility paths that must remain. The successful toy run used a
900-second cap. Failed attempts were retained and never published as latest.

Future canonical calibration runs export a run-local `best_params.py` only
after a fresh native check, exact repeat and matching complete input/grid
fingerprint. The export preserves effective fixed primitives, the verified price
and derived H0. Synthetic tests cover altered caller inputs and reject mismatched
or missing identity; no new optimization was run to test export. The historical
preparation's zero unsecured-credit binding is now explicit, exactly as in the
authenticated predecessor. Existing live cluster bundles were not changed.
