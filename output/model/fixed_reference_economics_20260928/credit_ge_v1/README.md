# Closed stationary borrowing equilibrium

**2007 stationary reference — block0506, September 28 verified export**

The author requests GE, with house prices and population endogenous. The
completed sibling `credit_v1` is a fixed-price mechanism, not this result.
Current state: bounded GE job **18764943** submitted after lead review and Torch
zero-solve syntax/closure/bracketing checks. No GE result certified yet. Sources
are immutable in `/scratch/td2248/projects/fixed_reference_credit_ge_20260929/sources_solve_v1/`;
local `source.sha256` pins them. Results and each completed-case receipt are in
the sibling `ge_results/solve_v1/`, with `solve_v1.log` at the remote root.
The q0 controller smoke passed: saved credit arrays, 14 fits and 31 parameter
estimates replay exactly, with all 17 standard plots. The price search has
started. Its 0.022947 birth-renewal residual is explicitly not an equilibrium;
receipt: `q0_smoke_receipt.json`.
Torch zero-solve
preflight 18763586 passed in 27.9 seconds at price factors 1, 1.05 and 1.35:
independent solvency recurrence, exact 262-node q0 grid, zero infeasible entrant
or inherited mass. Receipt: `preflight_ge_receipt.json`.

Remove the same artificial renter, purchaser and incumbent-owner debt limits,
retaining lifetime repayment and nonnegative net estates. Keep preferences
(including child benefit), earnings, interest, survival, taxes, housing menus,
entry endowments and the absolute housing-supply curve fixed. No recalibration.

Let B be adjusted births, E entrant households and d housing demand, all per
normalized household. Price q solves B/(2.1 E)=1. Population is then
N=Hs(q)/d(q), with Hs the unchanged absolute housing supply schedule. Verify
PAYGO at the actual stationary age/income distribution. The pension can remain
at its reference value only if that verification passes. Birth renewal, fiscal
balance and housing clearing must all hold; changing population alone is not GE.

Thus completed fertility returns to replacement through prices and population,
without changing child benefit. The endpoint comparison concerns population,
prices, rents, ownership and fertility timing. A transition and its temporary
fertility response remain separate work.

| Closure object | Classification and retained treatment |
|---|---|
| Preferences / supply intercept | Reference estimated or calibration-normalized values, held fixed |
| Earnings and survival | Externally estimated/fixed reference processes |
| Initial population / entrant endowments | Empirically normalized; preserve every original entrant atom |
| Birth-to-household conversion / entry clock | Author-fixed 1/2.1 conversion; half after 16 years, half after 20 |
| Geography | One pooled market, closed population; no outside entry |
| Housing supply | Fixed reference absolute curve, elasticity 0.63 |
| Payroll tax / pension | Reference tax fixed; actual stationary PAYGO must pass |
| Estates | Provisional reference net-estate funding and residual sink retained; counterparties/physical settlement remain outstanding |
| Price and population | Endogenous outputs, not normalized after the shock |

Independent mathematical review confirms the natural limit at any stationary
q: renter human-wealth limit is unchanged; an owner subtracts net housing
liquidation value. This is not a dated-price transition limit. Candidate grids
retain original160 nodes and entry atoms, adding the same economically derived
boundary knots at each price. Failed support or a required node outside the
retained range stops the run.

Preparation used worker_fast (Terra medium) for the adapter and initial driver.
The adapter passed; the driver draft failed lead review before any launch
(incorrect fields/accounting and incomplete root/repeat logic). Its rejected
source is retained separately. A single GPT-6 Sol coding repair has a 20-minute
cap and no further delegation; it reuses the completed credit driver. Lead
reviews critical changes. Torch preflight: zero solves, three prices, five
minutes. Proposed solve cap: 12 new lifecycle solves including
reference-price replay and exact selected repeat, 300 seconds per case,
40 minutes total, one CPU and 16 GiB. No hidden retries or deadline extensions.
Source/plan/checkpoint identities are pinned before launch. A failed or
unbracketed result is retained and reported, not turned into a certified root.

Keep full fit/parameter tables and the unchanged 17 diagnostics accessible;
the author-facing deliverable is one compact GE comparison, without a PDF.
Large checkpoints and runtime source snapshots stay on Torch.

Overnight supervision: `check-frozen-reference-borrowing-ge` checks every ten
minutes. The author explicitly authorized investigating hiccups rather than
only notifying. At most two targeted cheaper-worker repair cycles (20 minutes
each) may fix diagnosed implementation errors in new immutable versions, with
Torch verification and lead review. Preserve failed evidence and original
remaining numerical budgets/deadlines; no economic changes, relaxed gates,
unchanged retries, or duplicate searches. Pause after completion or a concrete
unresolved blocker, and no later than September 29 at 09:00 New York time.
