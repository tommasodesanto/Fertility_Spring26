# Borrowing mechanism at fixed prices

**2007 stationary reference — block0506, September 28 verified export**

Question: does access to borrowing against future earnings change births and
their timing? Remove artificial unsecured-credit, purchase/down-payment and
incumbent-owner debt limits. Retain all preferences (including psi), earnings,
prices/rents, taxes/pension, entry endowments, survival and repayment. This
first comparison is a fixed-price mechanism, not an equilibrium or transition.

Under the frozen zero-subsistence, zero-transfer, one-market specification,
the natural renter saving floor satisfies
\(L_j=\max\{0\text{ if death is possible},\max_{z'}(L_{j+1}-y_{j+1,z'})/R\}\),
with terminal floor zero. For owners subtract net housing liquidation value.
The maximum is over all income outcomes with positive probability. A separate
minimum-over-tenures recurrence independently verifies this reduction.

Torch preflight18753555 (29 seconds, zero solves) verifies the recurrence and
12 support fixtures. Every inherited household is economically solvent, but
the original160-node grid restricts age30 renter debt to0.9535 rather than
the economic1.6788 and excludes1.18e-8 inherited mass. No household solve was
launched on that grid. Preserve the failed readiness receipt in `preflight_v1/`.
Preflight18753911 passes (30seconds, zero solves):262nodes, zero inherited or
entrant infeasibility, every original point mass preserved exactly. Its age30
renter debt limit1.67754 approximates the economic1.67884; maximum tightening
across ages/products is0.0016005. This verifies support, not policy convergence.

The subsequent budget, conditional on a passing preflight, is one exact
reference control, a baseline on the refined grid, and one same-price credit
solve on that same grid: one Torch worker,16GiB,10minutes per case,30minutes
total, no identical retries. Compare credit to the baseline on the same grid;
report the original-to-refined baseline difference separately. A failed control,
solvency, budget, mass, estate-funding or occupied-value gate stops the loop.
Sources and grid are pinned before launch. Numerical grid changes are
disclosed separately from the credit change; they are not grid convergence.

Job18754242 was submitted once for this three-case comparison. Plan SHA256:
`60d9fd49670c450942bf780a18daf91b90b57c5dcebc7231ccc1639df7657ac8`.
Frozen sources are `source_solve_v1/` under the Torch root below. Expected
computation is roughly8–15minutes before queue delay; hard cap30minutes.
Household reoptimization and the economic comparison remain unverified until
all three case receipts pass. Do not treat the preflight as a credit result.

Keep the17 standard plots and complete14/31 tables as supporting files.
Deliver a compact birth/timing/tenure comparison; no assembled PDF.
All large checkpoints stay on Torch under
`/scratch/td2248/projects/fixed_reference_credit_20260929/`.
