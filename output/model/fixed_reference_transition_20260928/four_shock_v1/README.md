# One- and four-shock estimation

**September 28 launch authorization:** Tommaso authorized launch once ready, then
read-only overnight monitoring for major failures. This supersedes the earlier
test-only instruction. The authorized September 29 retry, current checks,
budgets and launch status are in [`launch_v3/README.md`](launch_v3/README.md).
The failed first launch remains in `launch_v2/`; preparation evidence follows.

The active driver is `code/model/tools/run_e5f_preference_estimation.py`.
The September 28 author correction supersedes the earlier announced-path
preparation: **estimate four successive surprises**, or estimate one permanent
2007 shock to the final 2020–2023 fertility window. Each surprise is believed
permanent when it arrives. Subsequent shocks are absent from that forecast.
The earlier slides' imposed one-shock path is background, not an estimate for
this reference. Historical estimates are never automatically transferred.

`e5f_preference_shock_fit.py` fits each scalar shock using solved fertility
residuals. Failed equilibrium, terminal or horizon checks cannot be scored.
The driver solves the endpoint internally at each candidate preference,
adjusting its price for demographic renewal without normalizing the preference.
It solves price and pension paths jointly, verifies a fresh final fit, then
implements the accepted forecast's first period. Both birth-entry queues and
the complete household distribution pass unchanged into the next surprise.
The first-period replay uses that vintage's next price and value function,
preserving the original expected rent even when the next shock revises prices.

The four retained targets are arithmetic means of published NCHS annual TFRs:
2008–2011, 2012–2015, 2016–2019 and 2020–2023. Their original CSVs and complete
target/measurement contract are hashed. The model uses the retained sum of
age-specific four-year birth-flow/household-mass rates; this household analogue
is not a literal female-exposure TFR. Four sequential fits use all four rows;
the one-shock fit uses only the final row and reports the others as validation.
All rows, shock estimates/bounds and unchanged reference 14/31 tables are saved.

The initial distribution, both birth-entry queues, earnings, preferences other
than the explicit shock, and saved DUE credit are inherited unchanged from the
**2007 stationary reference — block0506, September 28 verified export**.
Physical stock is held at its actual reference quantity, not its supply-curve
intercept. Prices and period pensions jointly clear housing and PAYGO under
perfect foresight; rents follow the dated asset-pricing relation. The saved
elastic housing rule remains an explicit alternative. No credit experiment or
alternative economic specification is adopted. The historical fit is now
authorized after readiness; the preparation checks below held preferences fixed.

## Current estimator verification

The test-only harness is `check_estimation.py`, submitted with
`check_estimation.sh`. It cannot call the historical shock optimizer. Its native
checks hold the saved preference fixed: endogenous endpoint, six-date baseline
root, first-period replay, then the remaining five dates from the carried state.
The latter must reproduce the original forecast, distribution and both queues.
**PASS: 34 synthetic tests and the unchanged-preference native integration.**
Job **18759951** completed the endpoint, six-date root and first-period replay,
then failed a test assertion that compared relative period counters without
their one-period offset. Every economic row already agreed. That failure is
retained; no model source or numerical tolerance was changed to repair it.
Job **18761094** reran the 34 tests and only the carried five-date replay, using
the saved certified endpoint/forecast/checkpoint. All economic rows, final
household cells and both queues reproduce exactly (maximum differences zero).
This replay used one Bellman solve and nine exact-cache hits; 129.502 seconds
including setup, one CPU, less than 5 GiB recorded peak batch memory.

`estimation_tests/18761094/readiness.json` pins the tested estimator sources,
links the original evidence and owns both disabled estimator plans. The original
test had a 40-call ceiling; the targeted recheck had a 10-call ceiling. Both
held the saved preference fixed, with no historical fitting. Baseline cache
reuse is not a speed measurement for a shocked path.

## Earlier inner-engine verification

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

## Numerical safeguards and future execution

`e5f_four_shock_acceleration.py` supplies a two-block physical Jacobian and
uniform step scaling. `measure_jacobian` constructs a fresh approximation from
five explicitly budgeted mappings; old three-block derivatives are rejected.
No native derivative measurements were launched. The exact policy cache is
bounded at 2 GiB. Production mappings persist dated audits, latest/best
checkpoints and a heartbeat; final acceptance requires a fresh replay, both
entry queues and terminal distribution/population convergence. Standard
diagnostics are wired for first/middle/final dates. Longer-horizon comparison
and visual review remain separate requirements.

New draft plans remain `execution_enabled: false`. For the authorized launch:

- Choose increasing terminal horizons and finite time/evaluation budgets.
  Each forecast must satisfy terminal and horizon checks before its fertility
  enters the objective. No shocked-horizon adequacy is claimed from baseline tests.
- Pin a separately enabled plan after current-source readiness passes. The
  author's September 28 instruction already authorizes the launch.

The provisional estate settlement remains an economic outstanding item in the
canonical status. These checks certify the unchanged native mapping and code
wiring, not a shocked equilibrium, fitted history or horizon convergence.

Inspect a new estimator draft without importing the model or launching a fit:

```sh
python code/model/tools/run_e5f_preference_estimation.py --plan output/model/fixed_reference_transition_20260928/four_shock_v1/estimation_tests/18761094/four_successive_draft_plan.json
python code/model/tools/run_e5f_preference_estimation.py --plan output/model/fixed_reference_transition_20260928/four_shock_v1/estimation_tests/18761094/one_permanent_draft_plan.json
```

The older `runs/18754009/*_draft_plan.json` plans described supplied announced
paths and are superseded for this task. Keep them with their original evidence.
All imports, tests, hashing and numerical work run on Torch. `check_engine.sh`
and `check_engine.py` reproduce only the bounded preparation suite against the
staged `source_v2/` files; they cannot launch a preference shock. The remote stage
is `/scratch/td2248/projects/fixed_reference_transition_20260928/four_shock_v1`.
Do not overwrite an in-use source stage. Only compact receipts return locally.
