# One- and four-shock transition preparation

The shared engine is `code/model/tools/run_e5f_preference_transition.py`.
The preferred case has four preference levels in 2007, 2011, 2015 and 2019,
all known in 2007. The alternative is one permanent change in 2007, as in the
transition plotted in the September/JMP slides. The last level persists.
The historical estimates are not automatically transferred to block0506.

The initial distribution, both birth-entry queues, earnings, preferences other
than the explicit shock, and saved DUE credit are inherited unchanged from the
**2007 stationary reference — block0506, September 28 verified export**.
Physical stock is held at its actual reference quantity, not its supply-curve
intercept. Prices and period pensions jointly clear housing and PAYGO under
perfect foresight; rents follow the dated asset-pricing relation. The saved
elastic housing rule remains an explicit alternative. No credit experiment,
preference re-estimation, endpoint solve or shocked transition was run.

## Verification

Torch job **18754009 PASS**, 185.809 seconds, one CPU, two actual household
solves across the suite. `runs/18754009/readiness.json` pins the exact sources.
The 21 pure tests cover announcement timing, launch guards, fixed stock,
cache semantics, two-block Jacobian ordering/scaling, uniform Newton step,
root replay, derivative construction and horizon comparison.

| Native no-shock check | Household solves / cache hits | Largest housing gap | Largest PAYGO gap |
| --- | --- | --- | --- |
| One date, saved elastic rule | 1 / 1 | 1.760e-9 | 3.172e-13 |
| Six dates, fixed physical stock | 1 / 11 | 5.692e-8 | 7.133e-8 |

The one-date numerical rows exactly match the earlier uncached result. Six
dates cover both 16/20-year entry lags; adjusted/raw queues are retained. Their
maximum changes are 2.933e-8 and 4.885e-9, respectively, with population L1
change 1.222e-7, consistent with the small retained renewal discrepancy.
The terminal population, normalized distribution, both queues, price and rent
also pass explicit 1e-6 gates. A failed raw queue now marks the terminal result
`not_converged`; its rejection and nonfinite-tolerance guard are tested.
All household, purchase, estate-funding, probability, occupied-value and mass
gates pass. Inherited feasibility projection is zero. Exact cache reuse here
benefits from a stationary continuation; it is not a measured speedup for a
shocked path.

The source review's omitted-dependency concern was checked: `load_reference`
calls `e5f_evening_calibration_runtime.setup`, which verifies all **1,241**
sources in the frozen manifest, including the native perfect-foresight driver,
estate audit and DUE audit. The five reused numerical overlay modules are
byte-identical to their frozen import locations. The new engine/helper have
their own source pins. `runs/18754009/source_review.json` records this check.
The frozen project is mounted read-only. Earlier job 18753571 passed the first
20 tests and both native mappings; the final repeat verifies the added terminal
checks and corrected raw-queue status, with unchanged numerical rows.

## Ready code and outstanding inputs

`e5f_four_shock_acceleration.py` supplies a two-block physical Jacobian and
uniform step scaling. `measure_jacobian` constructs a fresh approximation from
five explicitly budgeted mappings; old three-block derivatives are rejected.
No native derivative measurements were launched. The exact policy cache is
bounded at 2 GiB. Production mappings persist dated audits, latest/best
checkpoints and a heartbeat; final acceptance requires a fresh replay, both
entry queues and terminal distribution/population convergence. Standard
diagnostics are wired for first/middle/final dates. Longer-horizon comparison
and visual review remain separate requirements.

Both saved plans remain `execution_enabled: false`. Before a future launch:

- Specify the new-reference shock levels and their provenance.
- Supply the matching verified stationary endpoint, complete fit/parameter
  tables, numerical settings, horizon and finite time/evaluation budgets.
- Receive the author's subsequent instruction to run. The current instruction
  authorizes preparation only. No separate approval system is invented here.

The provisional estate settlement remains an economic outstanding item in the
canonical status. These checks certify the unchanged native mapping and code
wiring, not a shocked equilibrium, fitted history or horizon convergence.

Inspect either draft without importing the model or launching a transition:

```sh
python code/model/tools/run_e5f_preference_transition.py --plan output/model/fixed_reference_transition_20260928/four_shock_v1/runs/18754009/one_permanent_draft_plan.json
python code/model/tools/run_e5f_preference_transition.py --plan output/model/fixed_reference_transition_20260928/four_shock_v1/runs/18754009/four_announced_draft_plan.json
```

All imports, tests, hashing and numerical work ran on Torch. `check_engine.sh`
and `check_engine.py` reproduce only the bounded preparation suite against the
staged `source_v2/` files; they cannot launch a preference shock. The remote stage
is `/scratch/td2248/projects/fixed_reference_transition_20260928/four_shock_v1`.
Do not overwrite an in-use source stage. Only compact receipts return locally.
