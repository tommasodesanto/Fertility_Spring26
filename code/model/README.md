# Model guide: what to edit, run, and inspect

## Start here

Your main folder is `code/model/`. The everyday workflow is:

**Choose parameters → solve the steady state → inspect saved results.**

| File | What you use it for |
|---|---|
| [`parameters/best_params.py`](parameters/best_params.py) | The current working calibration and fixed inputs. |
| [`parameters/toy_params.py`](parameters/toy_params.py) | **Your editable experiment:** change beta, financing, borrowing limits, or other supported inputs. |
| [`run_model.py`](run_model.py) | **The file you run** to compute a complete steady state. Its `PARAMETER_FILE` line selects best or toy. |
| [`plot_model_policies.py`](plot_model_policies.py) | Plot household decisions across wealth, age, and family states. Edit this file to customize policy figures. |
| [`plot_model_aggregates.py`](plot_model_aggregates.py) | Plot population averages and lifecycle profiles. |
| [`tools/start_model_explorer.command`](tools/start_model_explorer.command) | Start the browser explorer for saved results. It does not solve the model. |
| [`tools/start_model_playground.command`](tools/start_model_playground.command) | Open a Python console for quick experiments at a fixed price. This does not compute a complete equilibrium. |

For your first experiment:

1. In `run_model.py`, set `PARAMETER_FILE = "toy_params.py"`.
2. Edit the values inside `toy_params.py`, such as `beta_annual` or `phi`.
3. Run `run_model.py` to solve and save the experiment.
4. Run either plotter or the explorer launcher. With no command-line override, they follow the same `PARAMETER_FILE` choice.

You do **not** need to run calibration to try your own parameter values. The
[commands below](#edit-inputs-and-run-one-case) also show how to select toy inputs
with `--params toy_params.py` without changing the default in `run_model.py`.

Your latest complete results are kept separately:

- **Working calibration:** [`output/model/local_solution/latest/`](../../output/model/local_solution/latest/).
- **Toy experiment:** [`output/model/experiments/toy_params/latest/`](../../output/model/experiments/toy_params/latest/).

The isolated [joint birth-count experiment](experiments/birth_count_choice/README.md) retains the current parameters and has its own complete-fit readout and output root.

Each contains `SUMMARY.md`, full target and parameter tables, and
`standard_diagnostics/`, `policy_plots/`, and `aggregate_plots/`. Plotting reads
these saved results without another solve. Parameter files isolate inputs and
outputs; changing the shared solver equations requires a separate checkout for
code isolation.

## Calibration versus the steady-state solver

The **steady-state solver**, [`production/equilibrium.py`](production/equilibrium.py),
takes fixed parameters and computes equilibrium. At each trial price it solves
household choices backward and propagates the distribution forward. It adjusts
the price to satisfy birth renewal; with the housing-supply coefficient \(H_0\)
fixed, population adjusts to clear housing. `run_model.py` calls this solver for you.

The **calibration search**, [`experiments/purchase_timing_sandbox/calibrate.py`](experiments/purchase_timing_sandbox/calibrate.py),
repeatedly changes parameter values, calls that same solver, and compares model
moments with empirical targets. Its adapter is [`production/calibration.py`](production/calibration.py).
Calibration uses the population-one normalization and derives \(H_0\). Future
canonical runs export a verified, run-local `best_params.py`; they do not
automatically replace your working default.

## Detailed code map and current reference

**Read [`../../CALIBRATION_STATUS.md`](../../CALIBRATION_STATUS.md) first for the live specification, calibration, overnight jobs, and unresolved numerical questions.** This file is navigation, not a second status ledger. The October 3 stationary-GE deployment commit is `94c0a6e3`; its verified sources, receipts and limits are in [`../../output/model/production_deployment_20261003/README.md`](../../output/model/production_deployment_20261003/README.md). The prior contents of this README are preserved byte-for-byte in [`../../calibration_archive/model_frontend_20261003/README_code_model_before_map.md`](../../calibration_archive/model_frontend_20261003/README_code_model_before_map.md), SHA-256 `aea46dc8c37697791d4c7f989640c2d9a9863acdda21c0603aaefe699d382ccd`.

## Current stationary GE

- [`run_model.py`](run_model.py) is the editable local general-equilibrium entry point. It calls the canonical implementation in [`production/`](production/README.md): selected inputs, equilibrium and price root, native household and forward-distribution engine, reporting, and storage. The current runner uses the fixed-H0 closure and the adopted post-interest, soft-constraint chain-13 input snapshot.
- Each successful solve with the canonical `best_params.py` publishes a validated case under `output/model/local_solution/cases/` and updates `output/model/local_solution/latest`; toy runs use their separate experiment folder. The case includes full fit and parameter tables, saved arrays, 17 standard diagnostic plots, 8 policy plots, 7 aggregate plots, and explorer assets.
- [`plot_model_policies.py`](plot_model_policies.py) and [`plot_model_aggregates.py`](plot_model_aggregates.py) read `latest` without another solve. The [production guide](production/README.md) gives the authenticated Python command, cache behavior, supported inputs, and exact numerical limits.
- The default local GE uses `fixed_h0`: it holds the physical housing-supply coefficient fixed while solving the birth-renewal price root. Its reported implied H0 and population scale N are computed algebraically from those same solved policies and price, without a second solve. The `population_one` calibration normalization instead fixes household scale at one and derives H0. This accounting relies on the documented conditional scale-independence and fixed-payroll assumptions; it does not certify another equilibrium.

| File | Responsibility |
|---|---|
| [`production/inputs.py`](production/inputs.py) | Loads and validates the selected primitive and parameter snapshot. |
| [`production/parameter_files.py`](production/parameter_files.py) | Reads data-only parameter files, validates their inputs, and routes their output directory. |
| [`production/equilibrium.py`](production/equilibrium.py) | Runs the stationary renewal-price search and its acceptance gates. |
| [`production/native_phase_b.py`](production/native_phase_b.py) | Implements the bounded stationary closure and phase-level equilibrium checks. |
| [`production/native_price.py`](production/native_price.py) | Solves the native household and distribution problem at one price. |
| [`production/engine/household.py`](production/engine/household.py) | Backward lifecycle household choices and policies. |
| [`production/engine/distribution.py`](production/engine/distribution.py) | Forward distribution propagation and aggregate measurement. |
| [`production/calibration.py`](production/calibration.py) | Checks the pinned target/weight contract and adapts calibration to the GE solver. |
| [`production/reporting.py`](production/reporting.py) | Binds the authenticated empirical observers and reporting tables. |
| [`production/workflow.py`](production/workflow.py) / [`production/storage.py`](production/storage.py) | Validates a local run, writes its case artifacts, and publishes the `latest` pointer only after successful storage checks. |

## Edit inputs and run one case

The editable, complete input files are [`parameters/best_params.py`](parameters/best_params.py) and [`parameters/toy_params.py`](parameters/toy_params.py). The first is the adopted working chain-13 point under the retained wealth target; it is not a global optimum. The toy file is a separate copy with annual beta lower by 0.001. Both files use the same engine and grid. Edit the toy copy for a personal experiment; parameter files separate inputs and outputs, not source code.

Each file sets ten `PARAMETERS`, fixed economic `EXTERNAL_INPUTS`, supported advanced `NATIVE_OVERRIDES`, a starting price, closure, and time budget. For example, edit `PARAMETERS["beta_annual"]` to change discounting. In the retained soft purchase rule, `EXTERNAL_INPUTS["phi"] = [0.75] * 4` selects a uniform financed share; its complementary 25% is the nominal non-financed share, not a strict liquid-cash threshold. The four values must match. `unsecured_credit_limit` controls renter debt capacity. Derived earnings and pension fields should not be edited directly; change the source earnings or payroll inputs. The file format accepts dictionaries, comments, literals, arithmetic, and list repetition, then validates fields and shapes before a solve.

From the project root, one command runs one stationary GE using the selected file; it does not search over parameters. Plotting and explorer commands inspect that file's saved result without solving again:

```sh
PROJECT=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
PYTHON="$PROJECT/output/model/publication_refactor_20260929/local_env_v1/venv313/bin/python"
"$PYTHON" "$PROJECT/code/model/run_model.py" --params toy_params.py
```

Pass the same `--params toy_params.py` option to the policy plotter, aggregate plotter, or explorer to inspect that experiment's saved case. Relative names resolve inside `code/model/parameters/`, regardless of the shell's current directory. Only the exact canonical `best_params.py` routes to `output/model/local_solution/`; most other files route to `output/model/experiments/<file-stem>/`. Noncanonical files named `best_params.py` receive a source-specific experiment directory to avoid collisions. A plotter or explorer with no cache for its selected file reports the missing case instead of opening another one. See [`parameters/README.md`](parameters/README.md) for the full input contract and [`production/README.md`](production/README.md) for solver and calibration details.

## Quick household inspection and explorer

For fixed-price partial-equilibrium experiments, launch [`start_model_playground.command`](tools/start_model_playground.command), which opens [`model_playground.py`](tools/model_playground.py) in the authenticated Python 3.13 environment. Its initialization loads the selected input file and matching saved case but does not solve. To choose an alternate file, use `start_model_playground.command --params toy_params.py`. If that file has no cached solution, the console reports that fact, leaves `sol=None`, and still creates `model` from the selected inputs for an explicit fixed-price solve. The [playground guide](tools/MODEL_PLAYGROUND.md) has the examples and interpretation limits.

```python
model.show_parameters()
model.params["beta_annual"] = 0.98
result = model.solve(price=0.72)  # one fixed-price lifecycle solve; no GE root
result.aggregates()
result.plot_policy(variable="consumption", age=30, income=4)
```

`model.params` exposes the ten displayed calibration coordinates. Ordinary supported primitive edits through `model.P` (or the documented input dictionaries) are carried into a solve. Parameter-mapped and recomputed fields should be changed through their source controls. Entry-distribution and structural grid edits are unsupported by this saved input snapshot and fail explicitly. A fixed-price solve does not clear housing markets, fit targets, or establish a new equilibrium.

The saved-case browser requires its local server. Run [`start_model_explorer.command`](tools/start_model_explorer.command), check that its configuration names the case you intend to inspect, and open the local URL it reports. It serves saved artifacts and performs no model solve. Opening the HTML file directly is not equivalent to starting the server; do not reuse an already-running explorer until confirming that its case matches the intended cache.

## Calibration, transitions, and experiments

- The isolated [normalized CES-limit share experiment](experiments/ces_normalized_shares/README.md) and [four-chain Torch workflow](../cluster/ces_normalized_shares_calibration/README.md) test the author-authorized childless share `alpha(0)=.733` and parent share `clip(.733-delta_alpha_jump-delta_alpha*m,.05,.95)`, with both parameters free on `[0,.25]`. The normalized denominator applies in all family states and `h_P=0`; the 11-free/11-scored experiment adds `family_rooms` (target 0.38509964969278165; weight 280.52808370152104). Old wealth target 6.92658379107299 and other economic inputs are retained. The starting birth cost is 1.9, still free on `[0,8]`. The final v5 three-context preflight and smoke 19132940 passed. Four six-hour chains are running as array 19133352; verified October 3 at 22:32:52 New York. See the [deployment packet](../../output/model/experiments/ces_normalized_shares/overnight_v1/README.md) for current receipts. No overnight calibration result or adoption is recorded.
- `production/calibration.py` validates the pinned target and weight fingerprints before solving and uses the selected chain-13 packet. The full fit has 14 rows and the parameter table has 31; see the live status and linked collection receipts for values and restrictions.
- A parameter-file run evaluates one input vector. The outer calibration search is [`experiments/purchase_timing_sandbox/calibrate.py`](experiments/purchase_timing_sandbox/calibrate.py): bounded Nelder–Mead changes the ten coordinates and calls the adapter in `production/calibration.py`, which evaluates each point with the same stationary-GE price-root solver under the calibration's `population_one` closure. Target scoring follows each GE solve; it is not the price-root criterion. A selected point is exported to a run-local `best_params.py` only after fresh native acceptance and an exact repeat. It does not replace the canonical adopted file automatically.
- Maintained dated-transition entry points are [`tools/e5f_current_transition_runtime.py`](tools/e5f_current_transition_runtime.py), [`tools/run_e5f_current_transition_smoke.py`](tools/run_e5f_current_transition_smoke.py), and [`tools/run_e5f_current_transition.py`](tools/run_e5f_current_transition.py). Their presence or importability does not certify a transition. Read the current status and the named transition packet before proposing or launching a run.
- Separate post-interest continuation searches use the retained old wealth target and an experimental new wealth target. They are distinct contracts; neither changes the working anchor by itself. Read [`CALIBRATION_STATUS.md`](../../CALIBRATION_STATUS.md) and its submission receipts for current array identities and queue state.
- Active and historical experiment code is indexed by its owning packet in `CALIBRATION_STATUS.md`; do not infer that an old README paragraph, script, or output folder represents current authorization or live state. Detailed development chronology formerly in this file is retained in the archived copy linked above.

## Retained research strands and archived code

The stationary `production/` package is the canonical local GE path; other
packages remain distinct research models, compatibility code, or diagnostics.
The older discrete-time model remains in [`dt_cp_model/`](dt_cp_model/), with
its May benchmark collection in the [archived benchmark index](../../calibration_archive/model_legacy_20261003/benchmarks/README.md).
The `intergen_*` and [`intergen_seq_fertility/`](intergen_seq_fertility/)
packages retain separate intergenerational/sequential implementations, while
[`experiments/`](experiments/), [`refactor_lab/`](refactor_lab/),
[`demographic_transition/`](demographic_transition/), and the maintained
transition entry points above own scoped prototypes, experiments, or dated
transition code. Timing-specific experiments remain under
[`experiments/purchase_timing_sandbox/`](experiments/purchase_timing_sandbox/).
None is an alternate default for `run_model.py`.

The Howard test and surrogate calibration code are archived with their READMEs
at [`intergen_housing_fertility_howard_test/`](../../calibration_archive/model_legacy_20261003/intergen_housing_fertility_howard_test/README.md)
and [`intergen_surrogate_calibration/`](../../calibration_archive/model_legacy_20261003/intergen_surrogate_calibration/README.md).
Frozen source authentication still pins ten Howard files, eight surrogate
files, and the original June runner. The old
[`code/model/intergen_housing_fertility_howard_test/`](intergen_housing_fertility_howard_test/)
and [`code/model/intergen_surrogate_calibration/`](intergen_surrogate_calibration/)
paths are compatibility symlinks into the archive; the original
[`code/model/run_intergen_model.py`](run_intergen_model.py) remains a regular
file because the authenticated historical runner depends on its source path.
Keep these references intact for source authentication and historical use;
none is a canonical stationary production engine. The pinned source list is
[`source_manifest.json`](../../output/model/overnight_calibration_20260928/contract_v1/source_manifest.json).
The old May [`PLAN.md`](../../calibration_archive/model_legacy_20261003/PLAN.md)
and the archived [`run_intergen_model.py`](../../calibration_archive/model_legacy_20261003/run_intergen_model.py)
remain available for historical use. Frozen observer, oracle, and reporting
artifacts under `output/model/fixed_reference_economics_20260928/` remain
read-only dependencies where current code imports or authenticates them.

## Copy-paste handoff

> Read `memory/AGENT_MEMORY.md`, the latest `memory/daily/` note, `CALIBRATION_STATUS.md`, and this model map before model work. Use the linked production, playground, or experiment packet for the task. Recheck mutable queue or run status from its live receipt; do not carry pasted runtime values forward as current facts.
