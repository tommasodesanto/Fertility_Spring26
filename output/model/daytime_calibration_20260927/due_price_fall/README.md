# Permanent 10% price-fall diagnostic

Both dated tests passed on September27. This is NOT an equilibrium transition.
Reference: original de_0093; all31 parameters fixed, no recalibration. Both arms
face a permanent10% lower house price and proportional user-cost rent, unchanged
fixed tax/pension, earnings, entry distribution, grid and exact originalg0.
DUE alone changes existing-owner debt treatment; purchase rules unchanged.
No default, insurance, redistribution or population reset.

Both pass original-state support, budget, purchase/stayer debt, dated estate
funding and population accounting. Negative gross/net estates and negative-estate
death mass are exactlyzero. Current estates fund actual next entrants from the
inherited split16/20 entry queue. Mass residual magnitude8.465e-16.
Household borrowing checks have zero violating mass; maximum stayer shortfall
8.882e-16 is floating-point dust. Fiscal residual approximately1.5e-13.

Dated raw birth flow0.120639 baseline versus0.120667 DUE (+0.023%). These are
flows on the same inherited distribution, NOT completed fertility or stationary
calibration moments. Housing residual10.316%/10.412% at imposed prices is reported
without claiming clearing. Full equilibrium and terminal/horizon remain unsolved.

Full dated comparison: dated_comparison.csv. All31 fixed parameters: parameters.csv.
No14-target calibration table is produced for a dated distribution. The prior
same-price stationary14-row comparison remains ../due_stayer_check/target_comparison.csv.
Both arms used one Bellman call, one at a time, single numerical threads; wall
time18.473s and30.533s inclusive of setup/export. No stationary distribution solve.
New immutable plan and parent process receipts are in the isolated worktree
output/model/daytime_calibration_20260927/due_price_fall/preparation_v1;
source/case paths and exact original-state hash are in summary.json.

Reviewed DUE core/accounting/date routing now integrated into main as explicit
optional functionality. Existing production contracts/defaults unchanged; a new
calibration source contract must explicitly activate DUE and origin-aware audits.
Supplemental stayer policies inspected in supplement/owner_conditional_policies.png.
6.715% of initial owners lie below the new collateral threshold;3.553% of actual
DUE stayers choose debt above it. Young indebted owners preserve consumption;
selected older policy curves mostly overlap. No realized stayer binds the
death-solvency floor in this experiment. These are selected conditional slices,
not a full policy certification. Standard17 stationary templates are not reused
with a dated payload.
