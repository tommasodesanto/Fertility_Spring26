# Run and plot the model in Python

Use three editable Python files for the usual workflow: [run_model.py](../run_model.py)
runs one local fixed-price stationary solve, [plot_model_policies.py](../plot_model_policies.py)
plots raw conditional policies from a saved run, and
[plot_model_aggregates.py](../plot_model_aggregates.py) plots population
aggregates from that same saved run. Run model first, then plot its saved result.
Open each file in the editor and press Play
with `code/model/.venv/bin/python` selected. The scripts resolve project paths
from their own locations, so they work from any current working directory.

The [verified workflow check](../../../output/model/fixed_reference_economics_20260928/model_control_scripts_v1/verification.json)
replayed all 11 checked baseline arrays exactly. The household solve took 6.2
seconds and the complete run, including standard figures, saving and validation,
took 9.9 seconds. Both plotters then loaded the saved result without solving and
wrote all eight policy and seven aggregate figures. These are measured times
for the default 120-node fixed-price run.

The equivalent terminal command uses the project interpreter and the absolute
project path:

```sh
PROJECT=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
"$PROJECT/code/model/.venv/bin/python" "$PROJECT/code/model/run_model.py"
```

The run driver starts at the selected soft experimental point, chain 16 / case
0046. Its default ten internal values are written explicitly in
`INTERNAL_PARAMETERS`; common native inputs are in `EXTERNAL_INPUTS`, and extra
native `P` edits go in `NATIVE_OVERRIDES`. The default fixed price is
`0.7266387868818555`. All four settings are editable at the top of
`run_model.py`, and the saved result records them together with the complete
effective native `P` object. Editing these settings creates a new experiment;
it does not update a calibration record or paper baseline. The selected point
is an experimental fixed-price starting point, not an adopted calibration.

Each run uses one CPU core and a 600-second time budget. It writes a unique
folder under `tmp/model_runs/<UTC timestamp>_<short id>/`, including
`native_result.npz`, `metadata.json`, and a `native_diagnostics/` packet with
the standard 17 PNGs and `summary.json`. The summary records the model's
completed-fertility measure and housing demand, supply, and residuals. The
fertility measure is descriptive; it is not an independent fertility renewal
gate. Run metadata records finite-value, probability-like-array, and
distribution checks, elapsed time, fixed price, and whether a general-equilibrium
price root or calibration was run. The completed-fertility statistic and
housing residuals are descriptive outputs from one fixed-price solve. The
driver does not solve the price root for market clearing and issues no
general-equilibrium convergence certificate.
After the saved result passes a full round-trip check,
`tmp/model_runs/latest.json` points to it. Each folder has its own identity and
remains available for later plots.

After a completed run, open and press Play on either plotter. Both load the
latest validated saved run and do not initialize reference parameters or solve
the model. To revisit an older run, set `RUN_DIRECTORY` at the top of the
plotter to that run's folder. If no completed run exists, first run
`run_model.py`.

The policy plotter writes eight PNGs under that run's `policy_plots/`:
consumption, next assets, renter housing, ownership probability, childless
first-birth attempt probability, consumption by income state, next assets by
income state, and consumption by family state. Defaults compare ages 26, 30,
and 34 at income states 3, 5, and 7 (counted from 1), use childless renter
policy and inherited-tenure branches, and show the central inherited-wealth
range. `WEALTH_LIMITS = None` shows the full grid; the default `"central"`
selects the 0.01–99.5 percent pooled-mass range plus one adjacent node at each
end. The
family-state chart separately includes `(n, m) = (0, 0), (1, 1), (2, 1),
(2, 2)`, where \(n\) is children ever born and \(m\) is children currently at
home. Settings also select owner size, buying or staying policy, asset axis,
wealth limits, and figure display. The curves use raw saved policies at original
asset-grid nodes; policy branches are not averaged or interpolated. Edit its
Matplotlib plotting calls to add a view that is useful for your question.

The aggregate plotter writes seven PNGs, `aggregates_by_age.csv`, and
`aggregates.json` under `aggregate_plots/`. It plots mean consumption, next
assets, inherited assets, financial-asset change, housing rooms, and ownership
by age, plus pooled inherited-asset mass. Consumption is a flow per model
period in mean annual gross-earnings units. Financial assets are stocks in the
same units. The asset-change series is \(b' - b\), with inherited \(b\) weighted
by the post-fertility, pre-tenure distribution; it includes housing transaction
cash flows and is not national-account saving. The wealth plot shows mass at
each saved asset node, not a density, so uneven grid spacing remains visible.
It defaults to the central 0.01–99.5 percent of the pooled mass plus one
adjacent node at either end; set `WEALTH_RANGE = "all"` to show the full grid.
Edit the explicit `ax.plot(ages, ...)` calls to customize an age profile. The
script writes ordinary Matplotlib figures, so you can add a horizontal line
with `ax.axhline(value, ...)` next to any selected plot.

The ten internal parameter values and their native meanings are:

| Parameter | Meaning |
|---|---|
| `beta_annual` | Annual discount factor; native period beta is `beta_annual ** period_years`. |
| `chi` | Owner housing-service premium. |
| `first_birth_fixed_cost` | Fixed utility cost of the first birth. |
| `kappa_fert` | First-birth choice shock/logit scale. |
| `kappa_fert_continuation` | Later-birth attempt choice shock/logit scale. |
| `theta0` | Bequest utility scale. |
| `h_P` | Physical room floor added at the first child. |
| `child_benefit_curvature` | Curvature of the child benefit by children at home. |
| `tenure_choice_kappa` | Tenure-choice logit scale. |
| `psi_child` | Child benefit scale. |

The native mapping sets `P.rho = P.rho_hat = 1/P.beta - 1` as discount-rate
fields; its default gross asset return is `P.R_gross = 1.08243216`. It
sets `P.eps_fert = P.kappa_fert`, `P.child_room_floor = True`,
`P.hbar_first_child_jump = h_P`, and `P.hbar_child_rooms = 0`. The other eight
values map to same-named `P` fields. Child utility uses
`psi_child * m ** (1 - child_benefit_curvature)`, where (m) is the number of
children currently at home. User edits are not clamped to calibration bounds;
the solver can reject invalid parameter combinations.

The fixed-price driver computes \(q = R_{gross} - 1\) and
`user_cost_rate = q + delta + tau_H` from its editable common inputs. It checks
that the retirement segment in `income` is constant before synchronizing the
pension fields. These derived inputs and the full effective `P` are saved with
the run.

## Optional interactive one-price work

The earlier interactive tool remains available at
`code/model/tools/start_model_playground.command`. It loads the authenticated
selected soft reference and initializes the solver, but it performs no solve
until `model.solve()` is called. Use it for direct Python inspection or a
one-price experiment:

```python
model.show_parameters()
model.params["beta_annual"] = 0.98
changed = model.solve()
changed.aggregates()
changed.plot_policy(variable="consumption", age=30, income=4)
fig, axes = changed.plot_aggregates()
fig, axes = changed.plot_aggregates(wealth_range="all")
```

`model.params` is the editable ten-parameter dictionary. `model.P` exposes the
initialized native parameter object for inspection and other direct edits; the
ten visible parameters take precedence if a field is changed in both places.
Edits to `P` are copied into each solve. `model.b_grid` and `model.solver`
expose the native grid and solver. Each `model.solve()` call solves the
household model once at the current fixed price; it does not recalibrate, solve
the renewal-price root, or clear the housing market. Set a fixed price
explicitly with `changed = model.solve(price=0.72)`.

To compare an experiment with the saved selected soft solution, without solving
the baseline again, run:

```python
baseline = model.saved_result()
changed = model.solve(overrides={"beta_annual": 0.98})
comparison = model.compare(baseline, changed)
comparison["overall"]
comparison["by_age"]
```

Use `model.reset()` to restore the initialized selected parameters and native
inputs after an experiment. Result objects keep their own parameter and
solution snapshots, so later edits do not change earlier comparisons.
`model_experiment.py` provides the same comparison as an editable script; its
`OVERRIDES` and `PRICE` settings are at the top.

For inspection, `changed.g` is the realized post-tenure cross-section,
`changed.g_stay_distribution` is its owner-stayer subset, and
`changed.g_beginning_distribution` is post-fertility and pre-tenure mass.
`changed.c_pol` and `changed.bp_pol` are conditional policies; owner stayers use
`changed.c_pol_stay` and `changed.bp_pol_stay`. Aggregates report consumption
per model period, financial assets in mean annual gross-earnings units, rooms,
ownership, the inherited asset distribution, and age profiles. The stationary
distribution follows the model's demographic construction; after a parameter
change it is not an equilibrium certification.

The selected starting point and price come from the authenticated soft
selection record at
`output/model/fixed_reference_economics_20260928/soft_timing_review_v1/soft_selected.json`.
The target and weight contract is checked against
`output/model/fixed_reference_economics_20260928/normalized_calibration_v2/`.
The active solver source is the pinned `small_credit_lab` overlay; its stages
are documented in `code/model/refactor_lab/README.md`. The child-benefit mapping
is in `code/model/refactor_lab/engine/child_preferences.py`.

The saved wealth grid has not been certified as converged. The focused
[asset-grid diagnosis](../../../output/model/fixed_reference_economics_20260928/asset_grid_diagnosis_v1/README.md)
finds measurable resolution error in some aggregates, while not establishing
full convergence for policies, targets, or equilibrium prices.
