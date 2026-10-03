# Python model code map

**Read [`../../CALIBRATION_STATUS.md`](../../CALIBRATION_STATUS.md) first for the live specification, calibration, overnight jobs, and unresolved numerical questions.** This file is navigation, not a second status ledger. The October 3 stationary-GE deployment commit is `94c0a6e3`; its verified sources, receipts and limits are in [`../../output/model/production_deployment_20261003/README.md`](../../output/model/production_deployment_20261003/README.md). The prior contents of this README are preserved byte-for-byte in [`../../calibration_archive/model_frontend_20261003/README_code_model_before_map.md`](../../calibration_archive/model_frontend_20261003/README_code_model_before_map.md), SHA-256 `aea46dc8c37697791d4c7f989640c2d9a9863acdda21c0603aaefe699d382ccd`.

## Current stationary GE

- [`run_model.py`](run_model.py) is the editable local general-equilibrium entry point. It calls the canonical implementation in [`production/`](production/README.md): selected inputs, equilibrium and price root, native household and forward-distribution engine, reporting, and storage. The current runner uses the fixed-H0 closure and the adopted post-interest, soft-constraint chain-13 input snapshot.
- Each successful solve publishes a validated case under `output/model/local_solution/cases/` and updates `output/model/local_solution/latest`. The case includes full fit and parameter tables, saved arrays, 17 standard diagnostic plots, 8 policy plots, 7 aggregate plots, and explorer assets. `latest` is the saved October 3 case until a later successful run publishes another.
- [`plot_model_policies.py`](plot_model_policies.py) and [`plot_model_aggregates.py`](plot_model_aggregates.py) read `latest` without another solve. The [production guide](production/README.md) gives the authenticated Python command, cache behavior, supported inputs, and exact numerical limits.
- The default local GE uses `fixed_h0`: it holds the physical housing-supply coefficient fixed while solving the birth-renewal price root. Its reported implied H0 and population scale N are computed algebraically from those same solved policies and price, without a second solve. The `population_one` calibration normalization instead fixes household scale at one and derives H0. This accounting relies on the documented conditional scale-independence and fixed-payroll assumptions; it does not certify another equilibrium.

| File | Responsibility |
|---|---|
| [`production/inputs.py`](production/inputs.py) | Loads and validates the selected primitive and parameter snapshot. |
| [`production/equilibrium.py`](production/equilibrium.py) | Runs the stationary renewal-price search and its acceptance gates. |
| [`production/native_phase_b.py`](production/native_phase_b.py) | Implements the bounded stationary closure and phase-level equilibrium checks. |
| [`production/native_price.py`](production/native_price.py) | Solves the native household and distribution problem at one price. |
| [`production/engine/household.py`](production/engine/household.py) | Backward lifecycle household choices and policies. |
| [`production/engine/distribution.py`](production/engine/distribution.py) | Forward distribution propagation and aggregate measurement. |
| [`production/calibration.py`](production/calibration.py) | Checks the pinned target/weight contract and adapts calibration to the GE solver. |
| [`production/reporting.py`](production/reporting.py) | Binds the authenticated empirical observers and reporting tables. |
| [`production/workflow.py`](production/workflow.py) / [`production/storage.py`](production/storage.py) | Validates a local run, writes its case artifacts, and publishes the `latest` pointer only after successful storage checks. |

## Quick household inspection and explorer

For fixed-price partial-equilibrium experiments, launch [`start_model_playground.command`](tools/start_model_playground.command), which opens [`model_playground.py`](tools/model_playground.py) in the authenticated Python 3.13 environment. Its initialization loads the production inputs and available saved case but does not solve. The [playground guide](tools/MODEL_PLAYGROUND.md) has the examples and interpretation limits.

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

- `production/calibration.py` validates the pinned target and weight fingerprints before solving and uses the selected chain-13 packet. The full fit has 14 rows and the parameter table has 31; see the live status and linked collection receipts for values and restrictions.
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
The old May [`PLAN.md`](../../calibration_archive/model_legacy_20261003/PLAN.md)
and [`run_intergen_model.py`](../../calibration_archive/model_legacy_20261003/run_intergen_model.py)
implementation are historical; [`code/model/run_intergen_model.py`](run_intergen_model.py)
is only a compatibility forwarder. Frozen observer, oracle, and reporting
artifacts under `output/model/fixed_reference_economics_20260928/` remain
read-only dependencies where current code imports or authenticates them.

## Copy-paste handoff

> Read `memory/AGENT_MEMORY.md`, the latest `memory/daily/` note, `CALIBRATION_STATUS.md`, and this model map before model work. Use the linked production, playground, or experiment packet for the task. Recheck mutable queue or run status from its live receipt; do not carry pasted runtime values forward as current facts.
