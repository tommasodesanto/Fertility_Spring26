# Borrowing mechanism at fixed prices

**2007 stationary reference — block0506, September 28 verified export**

## Economic result

At fixed prices, relaxing credit brings first births forward and increases
ownership substantially; completed fertility rises more modestly. All
artificial renter, purchaser and incumbent-owner limits are removed together,
so this does not separately identify the down-payment channel.

| Outcome | Baseline, matched grid | Solvency-only credit | Change |
|---|---:|---:|---:|
| Impact births per 1,000 households per four-year period | 115.51 | 122.49 | +6.04% |
| Impact first births per 1,000 households per period | 50.24 | 56.07 | +11.59% |
| Completed fertility, recomputed cohort | 2.1008 | 2.1482 | +2.26% |
| Mean first-birth age, recomputed cohort birth flow | 25.927 | 25.334 | -0.593 years |
| Ownership on impact | 67.00% | 72.80% | +5.80 percentage points |
| Ownership, recomputed cohort | 67.00% | 74.27% | +7.27 percentage points |

First births account for 83.4% of the immediate birth increase. Households
renting before the shock account for 95.6% of the increase; households with
nonpositive inherited net financial wealth account for 67.6%. These are
occupied-state contributions, not averages over hypothetical policy states.
Among renters after the tenure decision, the share with negative end-of-period
financial assets rises from 13.35% to 41.78% on impact. This measures debt
positions, not new loans; the baseline permits rolling inherited renter debt.

The 2.26% completed-fertility increase and earlier first births support a timing
and liquidity mechanism. Housing services rise 2.70% on impact but only 0.37%
in the recomputed cohort, alongside the much larger ownership change. Prices
and rents, earnings, preferences including psi, fiscal inputs and entry
endowments are fixed. This is a conditional cohort calculation, **not a new
market-clearing or demographic steady state, and not a transition**.

Supporting evidence, without an assembled report:

- [Comparison](solve_v1/comparison.json) and [completion receipt](solve_v1/completed.json).
- Full target-fit tables: [baseline](solve_v1/grid_control/target_fit.csv), [credit](solve_v1/credit/target_fit.csv).
- Full parameter tables: [baseline](solve_v1/grid_control/parameters.csv), [credit](solve_v1/credit/parameters.csv).
- Original parameter classifications and external restrictions: [frozen reference](../../fertility_identification_20260928/resume_v1/selected_export/primary/parameters.csv).
- The unchanged 17 standard plots: [baseline](solve_v1/grid_control/standard_diagnostics/), [credit](solve_v1/credit/standard_diagnostics/).
- Supplemental occupied-state [birth decomposition](summary_v1/impact_birth_decomposition.csv) and [debt positions](summary_v1/next_saving_debt.csv).

## Contract and verification

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

Job 18754242 completed exactly three solves in 517.60 seconds. Plan SHA256:
`60d9fd49670c450942bf780a18daf91b90b57c5dcebc7231ccc1639df7657ac8`.
Frozen sources are `source_solve_v1/` under the Torch root below. Expected
computation was roughly 8–15 minutes before queue delay; hard cap 30 minutes.
All three receipts passed. The original control exactly reproduced 113 arrays,
14 fit rows, 31 parameters and 17 PNGs. Repayment, occupied value monotonicity,
budget, mass, transaction and estate-funding checks passed; no occupied
unreachable continuation or negative estate exposure was found. Only 0.0283%
of cohort mass sits at a numerically tightened solvency floor. The matched-grid
comparison is not a full grid-convergence certificate. The original-to-refined
baseline completed-fertility change is +0.000786, reported separately.

Job 18755159 extracted the supplemental decomposition in 40.59 seconds, with
zero lifecycle solves and two supplied-policy period evaluations. Its fixed
budget was one CPU, 16 GiB and five minutes. Sources are `source_summary_v1/`
on Torch. Both jobs completed once, within budget, with no retry. To reproduce
the standard diagnostic set, use the frozen runtime's
`run_e5f_independent_numerical_audit.standard_diagnostics(packet, new_output,
validate_production_young=False)` on a saved remote case packet. Do not write
over these retained outputs.

Keep the17 standard plots and complete14/31 tables as supporting files.
Deliver a compact birth/timing/tenure comparison; no assembled PDF.
All large checkpoints stay on Torch under
`/scratch/td2248/projects/fixed_reference_credit_20260929/`.
Generated runtime source snapshots also remain on Torch; the local packet
retains their preparation receipts rather than duplicate source trees.
The exporter retains the legacy filename `lifecycle_2023.csv`; these files
describe the named 2007 reference and its credit experiment, not a 2023 result.
