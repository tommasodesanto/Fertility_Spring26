# Evening failure diagnosis — read-only, 28 September

Evidence is the collected final checkpoint of the completed evening search, plus its three selected-case stationary ledgers. This pass ran no model, imported no model modules, changed no source, and did not retry any case. Ordinary SSH initially failed, then recovered. A compact reduction of all 360 stationary ledgers is now saved in remote_ledger_summary.jsonl, with stage_counts.csv and remote_findings.json. Sample logs are also retained. Failed-price residuals/brackets were not saved by the original wrapper; explanations requiring those diagnostics remain hypotheses.

## What failed

The 360 search proposals produced 221 successes, 77 strict housing-market rejections, 61 owned 900-second timeouts, and one excluded late completion. Six final repeats are excluded from this count.

**All 77 housing rejections were broad-coverage proposals, and every one failed the first housing-equilibrium solve at the initial child-benefit intercept.** Thus they failed before fertility normalization could form a bracket. Broad coverage had 84 proposals: 77 rejected and seven timed out, with no success. These proposals independently sample all ten coordinates across their bounds (log scale where declared); they are not isolated perturbations of one parameter. There is consequently no identified single-coordinate cause.

No completed non-broad proposal failed the housing gate. The remaining 54 timeouts were 34 joint-local proposals, six first-birth fixed-cost coordinates, four child-benefit-curvature coordinates, three housing-taste-jump coordinates, three continuation-fertility-scale coordinates, two housing-supply coordinates, one owner-housing-premium coordinate, and one tenure-scale coordinate. Full counts and individual ten-coordinate vectors are in counts.csv and cases.csv.

By lane: block had 84 successes / 27 rejections / 9 timeouts; primary 77 / 24 / 19; identity 60 / 26 / 33 plus one late completion. Lane differences cannot be attributed mechanically to weighting: weights affect the evolving proposal center and therefore the economic points visited, not the equations solved at a fixed point.

## Numerical interpretation

The strict gate requires both solver convergence flags, with relative housing-demand-minus-supply error strictly below the model tolerance (at most 2.5e-5). It does not diagnose why the bounded price algorithm failed. The wrapper raises a generic RuntimeError and the normalization ledger saves that exception and elapsed time, but not the failed solution's scalar-refinement residual/bracket diagnostics. The checkpoint therefore cannot distinguish absent equilibrium from missing bracket, insufficient refinement, or demand discontinuities.

Code reading identifies a specific bounded-search issue worth checking, **not a confirmed cause**. The direct solver expands in the direction inferred from excess-demand sign. If no bracket appears, it resets a symmetric fallback bracket around the initial price, but retains the expansions counter. If the directional attempt exhausted max_expand, the symmetric bracket receives no further expansions. A root outside the searched interval, or nonmonotone demand, could therefore fail without establishing nonexistence. See intergen_eqscale_seq_optimized/solver.py, refine_one_market_markov_income, beginning around line 1633. A sign bracket also does not ensure the residual tolerance is attainable when discrete choices make demand discontinuous.

Warm prices carry only certified price/slope across child-benefit trials inside an objective. They do not warm-start the first trial from another candidate. Fertility normalization starts at the pinned intercept, probes the pinned step, doubles the distance until bracketed, then uses at most fourteen guarded secant steps; the objective additionally has a 23-equilibrium-call limit. The 900-second parent deadline can censor any stage. Recovered ledgers establish that every timeout had already completed four to eight certified equilibria: six cases completed four, twelve completed five, thirty-three completed six, eight completed seven, and two completed eight. Fifty-four were caught during a subsequent GE call. The remaining seven had completed normalization within its unchanged tolerance, as their final log lines confirm, but were killed during subsequent work before the objective receipt completed. These seven spent 846–864 seconds inside the completed GE calls; their exact remaining audit/export substage is not logged. No timeout happened in the first GE call.

For scale, the final primary winner used three stationary solves / eight price evaluations / 398.4 seconds; the block winner three / eight / 448.5 seconds; identity six / fifteen / 813.3 seconds. These successful examples make a 900-second censoring explanation plausible but do not establish it for any failed case. No conclusion that a failed point lacks a steady state is supported.

## Bounded matched retry recommendation, for lead review only

The failed-case stationary ledgers and selected stdout tails have now been recovered; preserve them and the old failures. Add an external diagnostic wrapper that records the returned GE price, residual, scalar_refine bracket_found/expansions/iterations, and warm-price trace before the unchanged strict gate is applied. Do not alter the economic equations, grids, targets, tolerances, or normalization rule.

At most six new numerical cases:

1. Two fixed-intercept GE diagnostic arms for broad rejected initial_0041_block: exact original starting price and, separately, half that starting price, keeping its ten coordinates and initial child-benefit intercept identical. Both arms retain all ordinary equilibrium and fiscal gates and a 900-second cap. This probes dependence on the price search region; it is not a recalibrated objective or a certified global-root test. Record complete price evidence even on rejection. If remote traces identify a different more informative case first, document the substitution before launch.
2. Two full-objective arms for timed-out joint-local initial_0044_block: identical settings with 900- and 1800-second parent caps.
3. Two full-objective arms for timed-out fixed-cost-coordinate initial_0097_identity: identical settings with 900- and 1800-second parent caps.

Maximum allocated case time is 7200 worker-seconds, at most three single-thread workers, with an explicit enclosing wall-clock cap (e.g. 60 minutes including startup/export), and the author deadline as the hard outer limit. Stop after these six; no automatic expansion. This is a numerical diagnostic budget, not a new search. The longer-cap arms test censoring only; passing requires the original complete scientific gate set. Compare exact parameters and target definitions, normalization output, residuals, solve counts, and timing; retain results as separately identified retries. A different price solution would require root/economics review rather than silent promotion.

Files: summary.json contains machine-readable totals and successful-case timing; counts.csv groups by design and lane; cases.csv lists all 360 proposals and all ten proposed coordinates. These are failure diagnostics, not a replacement target-fit report.

## Evidence update and sharper retry selection

The initial_0044_block timeout completed six equilibria in789.2seconds, then started a seventh at psi0.12124004. Its completed-fertility sequence crossed2.1 between psi0.12281100 (2.11922) and0.10281100 (1.87186), then refined toward the target. initial_0097_identity likewise completed six equilibria in832.1seconds and began a seventh. Neither sample shows a housing-convergence failure.

For a more distinct six-case design, replace the second timeout pair above (initial_0097_identity) with initial_0186_primary: it already reached completed fertility2.09999934 after seven equilibria/861.8seconds, then timed out after normalization. Thus the six cases probe one housing-search rejection, one still-normalizing timeout, and one post-normalization timeout. Keep original point vectors and all economic/scientific settings, caps and preservation rules above. This replacement is a recommendation, not a launched or approved run.
