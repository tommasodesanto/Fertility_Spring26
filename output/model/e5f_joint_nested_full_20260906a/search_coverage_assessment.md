# Lead assessment of initial search coverage

The independent review is confirmed against the actual 64 saved centers and
controller. All initial preference changes lie between -0.5091954181 and
-0.1577550125. At outer taste scale at most 0.1 they lie between -0.4364107182
and -0.2231046564. The initial design leaves near-zero declines unexplored.
Replacing rejected points with valid survivors does not retain this missing
joint region; differential evolution can extrapolate back but need not do so
within tonight's budget. This is a search-coverage issue, not evidence that any
parameter region or empirical target is infeasible.

The relevant within-tenure taste scale is lambda times kappa. Thus the concern
also applies to kappa=1, lambda=0.02, the actual rejected replay; a focus on low
outer kappa alone is too narrow. No quantitative derivative or required
preference decline has been established. A smaller decline is a proposal
heuristic to test, not an additional identifying restriction.

A possible unchanged-budget amendment would retain the exact anchor and all
47 original scale-grid cases, replacing the sixteen nearby random proposals
with sixteen paired scale-grid proposals. For each of eight original kappa
values and lambda in {0.02,0.2}, retain that original grid row's other nine
coordinates except the preference change, and multiply its physical preference
change by min(0.5,lambda*kappa/1.6). The denominator is the smoke anchor's inner
taste scale (2 times 0.8). This covers small changes without removing the
original large-change points or adding solves. The multiplier is capped at one half so every paired proposal differs from
its source. All eleven coordinates remain free in later
search and all bounds, targets, weights and gates stay unchanged.

This amendment is implemented LOCALLY AS PREPARATION ONLY, not yet adopted
for a cluster search. Twelve controller tests pass. All48 original anchor/grid
vectors remain exact, all16 new proposals are distinct and within the original
bounds, and every other controller function is AST-identical. The new paired
declines range from -0.1877535624 to -0.00003064868496. The full physical table
is proposed_initial_population.csv; proof is search_coverage_preparation.json. First inspect small-shock canary17091265's outcome/runtime and complete
smoke17090362. If the canary exposes a new defect, diagnose it first. If it
establishes unacceptable runtime, reconcile the case/time budget explicitly;
do not silently relax the three-timeout health rule. Any adopted amendment
requires reviewed controller/test changes, exact regeneration of the full
physical-coordinate table, new controller/source manifests and immutable
snapshot/contract, and actual imported-smoke preflight. Scientific bundle and
solver must remain unchanged; completed smoke can still certify the unchanged
scientific loop, with the new scheduler tested separately. No user permission
is needed for an evidence-backed starting-point change within authorized scope.

The normalization replay itself is fully adjudicated in
support_repair_e/replay_assessment.json; no further model repair is indicated
by that controlled rejection. The baseline anchors independently reproduce
all targets, parameters,253 historical numbers and17 PNGs exactly, taking
1969.84 and1970.42seconds. Current f helper/source are unchanged beneath the
running canary. The coverage worker has finished and owns no files.

Recorded: 2026-09-07T04:31:10.780330+00:00
