# Fixed unsecured-credit contract v1

Reference: `2007 stationary reference — block0506, September 28 verified export`.
The before manifest pins `parameters.py` at `66f86697c2c58ca3864305bf13dd2be71a008905b2beb573f1a4ebafabef5464`, `solver.py` at `b637a655a9344b63f4461ee0fa4796c04bd98188477c4e6ace2c48ae0fc8aec1`, and `kernels.py` at `639c9a21797dbc9f2a0e9a891f283c115353c2edfcb89c959a7fe9f32b86ca27`. `prepare_overlay.py` refuses either a hash mismatch or an existing destination.

`P.unsecured_credit_limit=None` keeps the legacy renter rollover-and-age-taper path. An explicit finite scalar \(D\geq0\) instead imposes renter saving \(b'\geq-D\), independently of current unsecured debt, age, income, and the existing taper arrays. At terminal age or with positive current death probability, the separate non-negative-estate restriction gives \(b'\geq\max\{-D,0\}=0\). Owners, purchases, legacy debt arrays, and mortality/estate rules are unchanged. `native_solvency_credit` and the scalar contract fail fast when combined; legacy/factored Bellman routes also fail fast for an explicit scalar, so the supported route is `solve_bellman_full_markov_income`.

An explicit scalar, including \(D=0\), turns on the owner-to-renter gate: raw liquidation wealth \(b+(1-\psi)pH\) must be non-negative before interpolation or price-grid clipping. The deterministic and logit tenure kernels and their Python fallbacks apply the same gate. The native renter saving kernel receives a separate scalar floor; no owner input array is changed.

`overrides.json` records the planned case \(D=0\); no positive magnitude is selected. The parent fixed-reference manifest SHA is `147f9e2cb20f66350f1ceaa16cb41f822041ec869676ef5d5b9d04f16e4190d4`; the checkpoint SHA `b15ba92dc60e3d5590d2beb6e05d36f71d17b20b1a432edc2c2db926a217309d` is recorded identity only and was not downloaded or rehashed. Preferences, entry, fiscal objects, and prices are inherited unchanged; there is no run-controller or GE acceptance.

From the repository root, run the pure checks with:

```sh
NUMBA_DISABLE_JIT=1 code/model/.venv/bin/python output/model/fixed_reference_economics_20260928/credit_no_taper_v1/fixed_credit_contract_v1/test_contract.py
```

They cover legacy `None` identity, zero and positive limits at several ages and
current balances, estate floors, invalid values, external overrides from old
parameter objects, direct rebuild validation, natural-credit conflict, sale
boundaries and an actual saving-kernel borrowing fixture. Owner kernels and call
sites pass AST comparisons. The separate lead checks below execute the
Python fallback and unsupported-route guard. No lifecycle, equilibrium or
checkpoint is read by these tests.

Known inherited-entry blocker: the earlier strict-\(D=0\) numerical run failed with two negative-wealth entrants. This packet deliberately does not change entry.

## Lead review and bounded compilation check

The lead reviewed every economic source change against the two constraints,
and independently executed the actual Python `_tenure_location_stage` fallback
from its source: deterministic and logit choices reject negative raw sale
balances and accept the exact zero boundary. The unsupported legacy Bellman
guard also executed successfully. This small check took 47 milliseconds and
performed zero lifecycle solves.

The first draft mistakenly combined the explicit floor with the old floor,
preventing positive borrowing from zero wealth. It was rejected before any model
run and is preserved in `overlay_initial_unaccepted/`. The corrected saving
kernel replaces the old renter floor when the scalar is supplied; its executable
fixture confirms a positive cap permits borrowing, while cap zero prohibits it.
Owner kernel and owner call sites pass exact AST comparisons against the reference.

Torch zero-lifecycle compiled smoke **18843355** was submitted once with a
five-minute limit, one CPU and 16 GiB. Source and test pins are checked before
execution. Remote sources are read-only; receipts go to `verification_v1/`.
This submission does not resume the expired numerical comparison budget.

The author subsequently permitted small local tests on the new laptop with one
core. Torch **18843355** was cancelled while pending (zero runtime), before the
local test started. The compiled local check **PASS** took 5.33 seconds wall time,
4.94 seconds CPU and 219,807,744 bytes peak resident memory. Numba, BLAS and OpenMP
thread limits were all one; wall and CPU limits were 180 seconds. Receipt and log:
[`verification_local_v1/receipt.json`](verification_local_v1/receipt.json),
[`test.log`](verification_local_v1/test.log). The tests assert that the actual
tenure and renter-saving dispatchers compiled. Frozen and effective source pins
and test hashes were checked before execution. Zero lifecycle solves and zero
checkpoint reads; full runtime integration, inherited-entry reconciliation and
a valid equilibrium remain uncomputed.

## Renewed overnight authorization, September 29

The author asked to continue and monitor overnight, retaining Torch for
substantial testing. This is new validation authorization; the earlier urgent
comparison deadline remains closed. A cheaper worker's partial
`runtime_validation_v1/` controller was not submitted and is retained as
incomplete. The reviewed `runtime_validation_v2/` smoke 18845938 failed before
Python because Slurm relocated the wrapper and its sibling path was relative
to `$0`; zero lifecycle evaluations occurred. Preserve v2. A narrowly corrected
`runtime_validation_v3/` uses explicit staged launcher paths: first a
zero-lifecycle checkpoint/import smoke, then at most one
fresh exact control with the overlay's scalar unset. One thread, 24 GiB,
300 seconds per control and 900 seconds total from its launcher entry; no retry.
No strict-zero lifecycle, price search or GE will be launched if checkpoint
affordability confirms the analytic entry blocker. The hourly heartbeat
`check-frozen-reference-borrowing-ge` now monitors this credit validation and
pauses after delivery/failure or by 09:00 New York September 30.

At the September 29 23:45 New York wake, v3 smoke **18846467 passed** in
41 seconds with zero lifecycle evaluations. Actual-checkpoint arithmetic
confirms the two strict-zero entrant failures, including buying infeasibility.
After lead receipt/source review, the one authorized scalar-unset exact control
was submitted as **18849552**. Do not duplicate it; pins, receipts and successor
paths are recorded in `runtime_validation_v3/README.md`. A valid revised-rule
GE still needs the author's entry/credit decision.

September 30 closeout: **18849552 passed** in 129 seconds, one lifecycle
evaluation; 113 native and 67 solution arrays, all 14 fit rows and all 31
parameter estimates match exactly. The 17 retained standard plot hashes
were verified. Full compact receipts/tables are in
`runtime_validation_v3/collected/control/`. The revised zero-credit GE remains
uncomputed due to inherited entrant infeasibility. The monitor is paused;
there are no pending numerical jobs or remaining runs in this packet.

## September 30 recommendation for author review — not adopted

Reference: **2007 stationary reference — block0506, September 28 verified export**.
The refactor handoff was assessed against its README, report and final
verification receipt; all 18 current engine hashes match that receipt.
`code/model/refactor_lab/` is suitable for the next isolated household and
distribution diagnostics. Completed certificates will not be rerun. Its
default equilibrium command clears a normalized household population against
an elastic supply curve; our closed borrowing experiment instead solves birth
renewal through price and scales population to housing supply. Reusing the
faster engine does not authorize switching that equilibrium closure, weakening
renewal gates, changing preferences or replacing the frozen namespace.

The entrant problem is economic, not just an asset-grid artifact. The PSID
entry input has two negative wealth/income-bin means, -2.22253 and -0.05264,
carrying 39.9803% of the model's primitive draw probability. This is not a
measurement of the survey's individual debt prevalence. The old builder
multiplies these ratios by annual gross income and interpolates onto the grid;
the current contract preserves that wealth marginal using a diagnostic income
rank coupling, explicitly not an estimated joint distribution. See
`code/data/psid_followup_mar2026/output/intergen_income_entry_targets_20260716/block2_entry_wealth_18_24.csv`,
`code/model/tools/e5f_earnings_wealth_contract.py:16`, and
`code/model/refactor_lab/engine/distribution.py:125`.

The grid includes exact zero; interpolation of a nonnegative entry point cannot
produce negative-node mass. The blocked node at -0.255814 represents a genuinely
negative underlying point near -0.177586. At that point, the two current-income
states have cash -0.049374 and 0.016907 before rent. The first remains
infeasible even at the exact point, so finer grids alone cannot remove the
problem. September 30 source correction: `c_min=0.04` and the 0.01 housing
floor are legacy output-only floors; the accepted exhaustive/exact-allocation
path reports its actual allocation without them (`engine/kernels.py:771`).
They are not economic minimum-spending requirements on that path. The second
off-grid point is therefore not proved infeasible by comparing cash with 0.04.
Only 0.00802196% of the entrant cohort is rejected by the frozen grid pairing;
that small figure does not justify zeroing all negative entrant wealth.

**Preferred next diagnostic if the author retains zero unsecured credit:**
test whether a minimally changed joint wealth–income pairing can preserve both
the frozen wealth and income marginals while assigning positive entrant mass
only to budget-feasible states. A mass-preserving transportation calculation
can test existence; native household feasibility must certify any candidate.
This would revise an unestimated correlation/selection assumption, not forgive
debt or remove households. It remains a proposal requiring the author's review;
no candidate pairing, economic run or parameter override has been adopted.
If it is impossible or empirically inappropriate, the entry contract must be
revised explicitly, or a positive constant unsecured limit must be justified
externally and selected by the author. A limit chosen solely to make the two
cells pass is not an empirical restriction.

The existing compiled tests establish implementation correctness, and the
scalar-unset replay establishes reference equivalence. Neither establishes a
corrected-credit equilibrium. No taper, sale-debt rollover, mass deletion,
entry truncation, transfer, positive credit amount, recalibration or fertility
normalization was introduced. The completed numerical budgets remain closed.

## Published borrowing-limit comparison — September 30

The published Boar–Gorea–Midrigan article and its official replication were
checked separately from earlier drafts. Their baseline uses a constant liquid
asset floor of -0.4 model units for renters and owners: 7.324% of average annual
income using the official README conversion. The bound is fixed outside the
estimated parameter vector. A separate empirical rationale for that value was
not found in the checked published text/appendix/replication documentation.
The 2017 draft's 3.6% value and tenth-percentile rationale must not be attributed
to the published calibration. [Source receipt and exact evidence](literature_published_v1/README.md).
This research does not adopt a new project credit limit or authorize a solve.
