# Lead review of the Fable audit

Fable completed on September 30 at 22:50 New York. Its report is an independent review, not a new calibration result. No recommendations have been applied to the running searches.

## Confirmed findings

- Candidate `psi_child` reaches the native parameters, derived coefficient and reporting checks. The parameter is estimated in the current experiment; the best reported point retaining its initial value does not mean it is fixed.
- The price start is hard-coded at 0.40056088974841536 in `../run_psi.py`. The executed price-search overlay reads its numerical caps separately from authenticated `q_ref`: [0.09873368771827608, 6.318956013969669]. A different first guess need not change those caps, preferences, entry or equilibrium equations. It can change the search path under a finite budget. The claimed twofold speed improvement has not been benchmarked. A change must also preserve reproducibility of the selected-point postcheck.
- In local chain 0, cases 0000 and 0007 differ only in the floor, from 1.02716049 to 1.12716049 rooms. Price falls from 0.67061721 to 0.63775276 (4.9006%); mean rooms rise by 0.21605886; the first-birth room response rises by 0.04086221. The renewal-price adjustment therefore materially couples fertility preferences and housing outcomes. This is evidence about the mechanism, not proof that the equilibrium closure is erroneous.
- The coefficient table's description still says “fixed benefit”; the actual coefficient uses candidate psi. That wording is stale metadata.

## Conclusions requiring qualification

- The assertion that credible calibration by morning is impossible is a prediction, not an established result. Limited evaluations and incomplete optimization warrant caution, not a declaration of failure.
- Small normalized kappa steps do not establish that those parameters are effectively unexplored. The code constructs physical steps equal to 10% of their seed values, then converts them to normalized coordinates. Chain 1's steps are 0.01992824 and 0.04893365.
- Four selected points cannot establish global invariance of early fertility, an unreachable room-response target, or rank deficiency. These require local derivatives or additional evidence. A pooled regression across different parameter regions, especially with endogenous price added as a separate regressor, is not a substitute for the parameter Jacobian.
- Square roots of loss contributions are weighted residuals. Calling all of them empirical standard-error units is unjustified without checking each weight's provenance; some scales are diagnostic conventions.
- There are ten scored moments, three zero-weight validation rows and one normalization row. The report's “four validation rows” is inaccurate.
- The static floor calculation is not a bound on a dynamic model with discrete owner choices, tenure selection, continuation births and endogenous prices. Rental caps may affect ownership, but the report does not establish their causal contribution to the current misfit.
- The audit did not finish tracing the empirical first-birth window, bequest specification or local/Torch observer-source equality. Its report discloses these limitations. Those points are not certified by the audit.

The exploratory objective deliberately defers plots and selected-price verification until the selected postcheck, as authorized. The current postcheck contract and economic acceptance gates remain unchanged. Proposed changes to tolerances, parameter bounds, targets, closure or live searches have not been adopted.

## Additional checks completed at 23:00 New York

The two current initial simplexes each contain a center and ten passed one-coordinate probes. The lead independently recomputed their singular values from stored base residuals. Both finite-difference matrices have numerical rank ten at relative thresholds 1e-6 and 1e-8. After normalizing each Jacobian column to unit Euclidean length, condition numbers are approximately 2,622 and 6,565: full numerical rank, but substantial local ill-conditioning. This does not certify statistical identification or stability to smaller steps. See [the saved-data calculation and provenance](lead_checks/initial_simplex_rank.json). No model calls were made.

The separate weak-direction calculation identifies $h_P$ and $\psi_{child}$ columns with correlation $-0.960$ and $-0.967$ at the two centers. The weakest step-scaled direction leaves nearly all predicted residual movement on early fertility at each center. These local derivatives are not a global or statistical identification result; see [the diagnostic](lead_checks/weak_directions.json).

Across 253 passed evaluations in the 23:00 monitor snapshot, age-25 children ever born ranges from 0.45757 to 0.58063, below the 0.80953 target. It is not invariant across the evaluated parameter vectors; these evaluations do not establish reachability or impossibility. The local and Torch copies of chain 1 case 0009 have identical parameter vectors, frozen observer identities and six observer source hashes. These checks resolve the audit's local/Torch observer-identity question for this matched point.

Evidence: [timestamped monitor snapshot](../deployment/monitor_snapshot/snapshot_20261001T030022Z/compact.json) and its companion source-comparison receipt.

Complete audit: [final.md](final.md). Existing full target and parameter tables remain in [the monitoring packet](../deployment/monitor_snapshot/REPORT.md). Additional checks using saved evaluations are recorded separately; they do not run the model.
