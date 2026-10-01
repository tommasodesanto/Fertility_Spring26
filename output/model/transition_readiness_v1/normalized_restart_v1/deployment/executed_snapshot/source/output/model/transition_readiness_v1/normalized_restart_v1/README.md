# One-shock transition from the normalized calibration

## Verified facts

The author authorized restarting the one-permanent-shock transition on October 1.
The selected reference is normalized-calibration chain 20, case `0028_nm`:
initial physical population one, housing coefficient `H0=6.851575289344519`,
price `0.7167873404451099`, and child preference `psi=0.17198899419542374`.
Its original-weight calibration loss is `30.371887956158005`; all 14 target rows
and 31 parameter rows match the selected native repeat, and its 17 standard
plot hashes are verified. Numerical repeat verification does not establish
optimizer convergence or identification. This point is selected for the
transition experiment, not adopted as a production calibration.

Complete [target fit](native_reports/ROOT/target_fit.csv),
[parameters with bounds](native_reports/ROOT/parameters.csv),
[handoff and contract](handoff.json), and
[standard diagnostics](native_reports/REPEAT/standard_diagnostics/) are retained.
The first-child housing floor is `2.3`, at its upper search bound. The lead
independently checked complete table identities, loss summation, parameter
bounds and two plots (housing clearing and fertility by age).

## Economic changes and experiment

Relative to the cancelled `actualfit_v3` reference, the author-directed baseline
contract sets population to one and internally calibrates housing supply.
The selected experimental recalibration also updates `beta_annual`, `chi`,
`first_birth_fixed_cost`, `kappa_fert`, `kappa_fert_continuation`, `theta0`,
`child_benefit_curvature`, `tenure_choice_kappa`, and `psi_child`, together with
the derived child-benefit CRRA coefficient. All estimates, restrictions, bounds
and statuses are in the complete parameter table. Earnings, initial income and
wealth distributions, entry, timing, transfers, housing floors and empirical
target definitions retain the named calibration contract. These are baseline
recalibration changes; the transition experiment itself applies one permanent
2007 child-preference shock and estimates its level from the original 2020–2023
target `1.64575`. The earlier three windows remain zero-weight validation rows.

After calibration, `H0` is fixed throughout stationary endpoints and the path;
physical population may change. The calibration's housing-inversion observer
is never installed in transition callbacks. Fresh reconstruction of the
selected initial distribution and both entrant queues, followed by fresh
measured derivatives, is required. No cancelled-run root, Jacobian, shock
estimate or 2023 checkpoint is accepted as evidence for this baseline.

## Verification and budget

Eighteen focused zero-solve runtime tests pass independently. The normalized
handoff authenticates 223 executed sources, complete tables and all 17 plots.
The runtime binds authenticated `H0` once and enforces the original `[0.2,80]`
restriction and population-one contract. Its comparison of calibration and
fixed-supply reference reports checks every target/parameter field and common
closure field within `1e-10`, allowing only explicitly declared supply-role
metadata changes; all 17 plot hashes must remain exact. The two fresh native
reference reports must also match each other exactly. Those actual numerical
reconstruction checks have not yet run at this preparation snapshot.

The author-authorized diagnostic comparison tolerances remain `0.01`, with
original strict comparisons retained. Fiscal, housing, stationary renewal,
replay, state and terminal gates are unchanged. Full-state, terminal and
production certification remain separate requirements. The hard deadline
remains October 1 at 20:24:47 New York. The conservative remaining cap is
1,343 native calls; see [prior-call accounting](prior_call_budget.json).
The deadline may bind before a fit. Submission and actual-start receipts belong
in `deployment/`; a preparation receipt does not establish a running job.
