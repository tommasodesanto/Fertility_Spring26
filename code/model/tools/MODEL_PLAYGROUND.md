# Python model playground

Start by double-clicking `code/model/tools/start_model_playground.command`. The
interactive session loads the saved solution and the authenticated selected soft
reference, but runs no model solve. It displays the ten selected primitives.
For one core, one explicit fixed-price experiment, use:

```python
model.show_parameters()
model.params["beta_annual"] = 0.98
changed = model.solve()
changed.aggregates()
changed.plot_policy(variable="consumption", age=30, income=4)
fig, axes = changed.plot_aggregates()
fig, axes = changed.plot_aggregates(wealth_range="all")
```

`model.params` is the editable parameter dictionary. `model.P` exposes the
initialized native parameter object for inspection and direct edits to other
native inputs; the ten visible parameters take precedence if the same fields
are changed in both places. Such `P` edits are copied into each solve.
`model.b_grid` and `model.solver` expose the native grid and solver. The ten
editable values and their meanings are:

| Parameter | Meaning |
|---|---|
| `beta_annual` | Annual discount factor; native `P.beta = beta_annual ** P.period_years`. |
| `chi` | Owner housing-service premium. |
| `first_birth_fixed_cost` | Fixed utility cost of the first birth. |
| `kappa_fert` | First-birth choice shock scale used in the logit. |
| `kappa_fert_continuation` | Choice shock scale for later-birth attempts. |
| `theta0` | Bequest utility scale. |
| `h_P` | Physical room floor added at the first child. |
| `child_benefit_curvature` | Curvature of the child benefit by children at home. |
| `tenure_choice_kappa` | Tenure-choice logit scale. |
| `psi_child` | Child benefit scale. |

The native mapping also sets `P.rho = P.rho_hat = 1/P.beta - 1` as discount-rate
fields; it leaves gross asset return `P.R_gross = 1.08243216` unchanged.
It sets `P.eps_fert = P.kappa_fert`, `P.child_room_floor = True`,
`P.hbar_first_child_jump = h_P`, and `P.hbar_child_rooms = 0`. The other eight
values map to same-named `P` fields. Every solve starts from a copy of current
`model.P`, applies the visible dictionary, and rebuilds the shared solver
objects. It starts from the authenticated selection unless you directly edit
other `P` fields. Values are not silently clamped to calibration search bounds;
the model can still reject an economically or numerically invalid value.
The child preference routine uses the selected curvature and scale as
`psi_child * m ** (1 - child_benefit_curvature)`, where `m` is children at home.

Each call to `model.solve()` runs the stationary household model once at the
current fixed price. It does not recalibrate, solve the renewal-price root, or
clear the housing market. Set a different fixed price explicitly with
`changed = model.solve(price=0.72)`. To compare against the saved selected soft
solution, which requires no baseline solve, run:

```python
baseline = model.saved_result()
changed = model.solve(overrides={"beta_annual": 0.98})
comparison = model.compare(baseline, changed)
comparison["overall"]
comparison["by_age"]
```

Use `model.reset()` to restore the initialized selected parameters and native
inputs after an experiment. Result objects keep their own parameter and solution
snapshots, so later edits do not change earlier comparisons.

`model_experiment.py` provides the same comparison as a readable script. Edit
its `OVERRIDES` and `PRICE` at the top, then run
`code/model/.venv/bin/python code/model/tools/model_experiment.py`. Importing the
script does not solve anything; executing it runs the changed case and displays
baseline and changed age profiles. The launcher does not invoke that experiment.

Use `changed.g`, `changed.g_stay_distribution`, and
`changed.g_beginning_distribution` for realized, stayer, and post-fertility
pre-tenure mass. `changed.c_pol` and `changed.bp_pol` are conditional policies;
owner stayers have `changed.c_pol_stay` and `changed.bp_pol_stay`. The aggregate
summary reports consumption per model period, assets in mean annual
gross-earnings units, rooms, and ownership. It also returns the inherited asset
distribution and age profiles. The stationary distribution uses the model's
demographic construction; it is not a calibrated equilibrium after a parameter
change. A changed parameter or price can alter the entire stationary
distribution, not just the plotted policy.
The aggregate plot defaults to the central inherited-asset range; pass
`wealth_range="all"` to show every saved asset-grid node.

The saved wealth grid has not been certified as converged. The focused
[asset-grid diagnosis](../../../output/model/fixed_reference_economics_20260928/asset_grid_diagnosis_v1/README.md)
finds measurable resolution error in some aggregates, while not establishing
full convergence for policies, targets, or equilibrium prices.

For one-off experiments, keep the workload to one core and a ten-minute budget.
The standard launcher sets Numba, BLAS, OpenMP, and NumExpr thread counts to one.
The saved wealth grid has 120 nodes from -12 to 3000. Zero mass at its upper
endpoint does not establish grid convergence.

The selected parameter record is authenticated before edits from
`output/model/fixed_reference_economics_20260928/soft_timing_review_v1/soft_selected.json`
and its pinned source checkpoint. The target and weight contract is checked
against `output/model/fixed_reference_economics_20260928/normalized_calibration_v2/`.
Saved solution arrays are loaded from the hash-checked cases named in
`output/model/fixed_reference_economics_20260928/soft_timing_review_v1/explorer_cases.json`.
The active solver is the pinned `small_credit_lab` source overlay; the model
stages are documented in `code/model/refactor_lab/README.md`. The child benefit
mapping is implemented in `code/model/refactor_lab/engine/child_preferences.py`.
