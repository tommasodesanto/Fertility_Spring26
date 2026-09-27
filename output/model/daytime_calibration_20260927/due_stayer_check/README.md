# DUE existing-owner rule: matched fixed-price check

The author selected the rule subject to a matched numerical check. The same-price
comparison now passes. The implementation remains isolated and default-off in
main until production integration; a dated price-fall check remains outstanding.
This is NOT a newly cleared equilibrium or a recalibration.

Reference: original overnight de_0093, not the later local-search candidate.
All31 economic parameters and prices fixed. The only economic change in the DUE
arm is grandfathering existing-owner debt: b_next >= min(b,-phi*p*H), with
separate current net-estate solvency when death is possible. Purchases retain
origination rules; interest is serviced. No PTI, transfer, default or insurance.

Complete14-row target/gap/weight/contribution comparison: target_comparison.csv.
All31 fixed parameters with inherited restrictions/bounds: parameters.csv.
Machine-readable checks and source paths: summary.json.

DUE household solve31.686seconds. Housing residual3.375e-4 at fixed reference
price, so the production market-clearing gate is NOT claimed. Budgets and debt
checks report zero violating mass; zero negative estates. Distribution nesting
L1=3.180e-15 and one-step L1=9.194e-14. Both17-plot standard packets inspected;
old-age ownership and wide-grid policy-tail features persist. No global policy
certification. Weighted diagnostic loss42.282 ->42.042, not optimized losses.

Reproducible source and full snapshots/plots/cases are under the isolated
worktree /Users/tommasodesanto/.codex/worktrees/due-existing-owner-credit/Fertility_Spring26,
output/model/daytime_calibration_20260927/due_stayer_check.
The bounded schedule retained its original absolute end1790536354.173106;
four model calls total, one worker at a time, single numerical threads. No
numerical worker remains and no price-fall arm was run.

Failures preserved: baseline_v1 reporting audit rejected an inherited tiny
infeasible tail; baseline_v2 weighted-loss absolute assertion was inappropriate
for weighted scaling, so an independent saved-state review verified exact core
arrays,14 model moments within1e-12 and actual weighted-loss formulas, then
exported17plots without another model solve. due_v1 exposed omitted stayer
saving in stationary reconstruction; due_v2 fixes this and passes.40 focused
worktree tests passed across the reviewed changes. No old failure was erased
or relabeled as a corrected run.

The initial-state guard separately forbids every redistribution while retaining
and explicitly recording the existing1e-12 numerical feasibility tolerance.
Saved-state replay preserves the baseline array exactly and rejects the larger
1.218e-8 frictionless transition support loss. That full transition is unsolved.

## Independent review follow-up, 15:29 EDT

Saved due_v2 all-origin negative estates are exactly zero. A missing explicit
rejection in the matched runner is now corrected for future source contracts,
using the production transition tolerance1e-10. Existing evidence is unchanged.
The standard17 plots show buyer-conditional saving/consumption and cannot alone
certify DUE stayer policies. A separate permanent10% price-fall dated comparison
is being prepared with exact originalg0; no new numerical run yet.
