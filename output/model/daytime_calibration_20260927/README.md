# September 27 daytime continuation and diagnostics

Author authorization: continue the calibration search; do policy, stationary
population, Jacobian and three-household-life diagnostics locally; compute a
same-parameter benchmark without the borrowing constraint. This supersedes the
completed overnight schedule's launch stop; the overnight heartbeat remains
paused. The reference is overnight continuation de_0093, with the exact source,
parameters and targets recorded in ../supervised_calibration_20260927/primary_final/.

Separate work and ownership:

- `search/`: Torch continuation, at most 24 workers and 192 search cases,
  two exact-loop smokes and two final repeats. Absolute 3.5-hour budget leaves
  90 minutes for verification/export. Economic model, targets, weights, bounds
  and grids unchanged; new starting point/seed and deadline are search changes.
- `jacobian/`: local central differences for all nine searched coordinates,
  18 objectives and at most two workers, 90-minute total/15-minute case caps.
  Each objective retains the endogenous child-benefit normalization. Weighted
  local sensitivity is not proof of global identification. Step definitions
  and fresh immutable fingerprints will be saved with the diagnostic.
- `households/`: saved-policy and full stationary-distribution checks,
  supplemental zoom plots and three reproducible illustrative household lives.
  No new equilibrium solution or change to the standard 17 diagnostics.
- Borrowing benchmark: separate experiment, never promoted into calibration.
  All preference parameters, including the normalized child-benefit level,
  must remain fixed. The author is clarifying whether to remove both the
  purchase and debt restrictions or just the down-payment restriction.
  Source mapping and a solvency-consistent definition precede the benchmark.
  Do not silently substitute a 100%-LTV experiment for unrestricted borrowing.

Torch authentication succeeds at the start of this session. Local diagnostics
have a shared cap of three simultaneous model evaluators (two Jacobian plus
one benchmark); saved-checkpoint work uses one additional process without
solves. Every new experiment writes to a fresh folder; overnight contracts,
checkpoints and final results remain immutable. All numerical failures are
preserved. Full fits and parameters accompany every new reported calibration.

Half-hour supervision is active on the existing task heartbeat, with these new
budgets. It will pause again after the authorized jobs finish and are reported.

Credit source review: the authenticated purchase-income adapter enforces
`x+y/R >= -phi*pH` and owner saving `b_next >= -phi*pH`; changing phi to one
relaxes both coherently but retains renter debt limits. Both purchase tests are
algebraically equivalent on supported states, so deleting just one is inert.
The estate audit explicitly leaves negative-estate creditor treatment unresolved.
A sure-repayment benchmark must respect liquidation solvency at possible death
dates; it cannot simply borrow to the numerical grid minimum. Full review was
read-only and did not change any economic mechanism.
