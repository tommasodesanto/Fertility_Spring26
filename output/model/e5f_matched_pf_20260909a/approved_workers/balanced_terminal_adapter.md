# Balanced terminal adapter: implementation handoff

September 11, 2026. Bounded implementation completed in the isolated checkout
`tmp/e5f_matched_pf`; no model solve, Numba compilation, cluster job, commit,
push, serialization or manuscript edit was performed by this worker.

## Owned deliverables

- `code/model/tools/e5f_balanced_terminal.py`
- `code/model/tools/test_e5f_balanced_terminal.py`

Both paths above are relative to the isolated checkout. This handoff is the
only worker-written file in the main checkout. Concurrent edits by other
workers were left untouched. The lead must review this adapter against the
mathematical specification before numerical use.

## Interface and economics

`solve_balanced_terminal(...)` requires the caller's parameters, grid, original
`supply_rule`, frozen `demographic_primitives`, `EndpointControls`,
`TerminalAuditControls`, initial asset price and **period-unit** pension,
price/pension bounds, separate fiscal tolerance, root slopes, log-step/damping,
maximum evaluations, monotonic deadline, conditioning/worsening limits, replay
tolerance and callback. The callback may explicitly be `None`. Only the optional
initial two-by-two Jacobian has a default. No population, migration, headship,
entry or housing-supply object is generated from diagnostic defaults.

The adapter requires the approved sequential architecture, a fixed payroll tax
of 0.179, four-year flows, period property tax 0.04 (1% annually), and zero
property-tax rebate. It preserves the caller's fertility preferences, utility,
other demographic primitives and supply elasticity/level. It does not rebase
supply, normalize fertility or prepare empirical inputs.

For every trial, including the reserved final replay:

1. Copy parameters and call `bind_social_security_income` with the trial
   period pension before household precomputation and the full
   `solve_markov_income_at_prices` call. No annual income resolver is used.
2. Reconstruct the stationary pre-fertility distribution as an inner seed.
3. Pass that seed, the trial policy, and the unchanged demographic and supply
   objects to `evaluate_endpoint`, which invokes the existing household/person
   stationary fixed point.
4. Re-evaluate the returned terminal pre-choice household distribution under
   that policy. Integrate Social Security over its actual post-choice household
   mass using `fiscal_accounts`; the seed age distribution and resident
   nonheads are excluded. Check that age/income mass is preserved and that the
   re-evaluated housing aggregates match the endpoint.
5. Send `(demand-supply)/supply` and
   `(payroll revenue-pension outlays)/max(revenue,outlays,1e-12)` to the existing
   `solve_social_security_path` with one date and fixed-tax closure. When both
   budget sides are zero, the residual is zero; when revenue is zero and outlays
   are positive, the residual is -1 and fails the gate. Property-tax revenue
   remains a separately reported unrebated surplus, outside PAYGO.

The endpoint mapping includes the existing population, one-step, birth-rate,
head, housing-law and finite-renewal gates. Household checks reuse
`primitive.dated_budget`, the existing numerical audit's occupied adjacent-wealth
value screen, and the matched-PF primitive probability bounds. Calendar
evaluation executes the core feasibility/dead-mass gates. All policy arrays are
checked for finiteness; reconstruction and projection masses retain explicit
limits. The adapter does not produce the 17-graph packet per trial.

## Memory and replay behavior

`BalancedTerminalResult.endpoint` contains the actual endpoint, policy and
parameters. `root_receipt` contains small scalars/dictionaries and root vectors,
with no policy or distribution tensors. At most the best and most recent
completed endpoint objects are retained outside the root. The generic root's
fresh final evaluation is counted within `max_evaluations`; the adapter returns
that exact newly evaluated object after matching its trial index, price and
pension. It does not recover a stale cached endpoint or add an unbudgeted solve.

`endpoint_production_eligible` requires generic root convergence, all mapping,
housing and Social Security gates, both residual replay gates, and association
with the actual fresh final endpoint. Failed or unfinished roots remain
diagnostic, including when the best prior trial passed. This flag establishes
only numerical terminal-endpoint eligibility: source/target contracts,
calibration, historical-path and horizon certification are still external
prerequisites for production use.

The deadline is checked between evaluations; an in-flight Bellman/person solve
requires the caller's process watchdog. This adapter owns no launcher, sizing,
wall-clock enforcement thread or output serialization.

## Verification

Pure local tests passed: **39 tests in 0.076 seconds**, including 12 new adapter
tests and the existing Social Security root and explicit endpoint test suites.
The test process used `NUMBA_DISABLE_JIT=1` and `PYTHONDONTWRITEBYTECODE=1`.
No economic runtime was loaded by these tests: the actual trial routing test
injects fake model/calendar objects.

Coverage includes fresh replay at an already cleared point, updates to both
unknowns, fiscal replay drift, invalid initial/final household mappings, bounded
evaluations, interrupted replay, fiscal/architecture/input rejection before
model work, zero-budget handling, exact income-binding order, single period
conversion, actual terminal versus seed accounting, unchanged supply and
demographic objects, and finite/probability/value/projection gates. A synthetic
tensor that raises on deepcopy verifies that root bookkeeping does not copy
endpoint tensors.

Reproduction from the isolated checkout:

```sh
NUMBA_DISABLE_JIT=1 PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=code/model/tools /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/.venv/bin/python -m unittest test_e5f_balanced_terminal test_e5f_social_security_root test_e5f_matched_pf_endpoint -q
```

No compiled fixed-price, actual terminal fixed-point, root convergence, runtime,
memory or production-grid result is claimed. The next step is the lead's source
review followed by the separately sized, immutable compiled exact-loop smoke.
