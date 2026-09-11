# Lead adjudication of the Fable 5.1 review

**Completed before any model implementation or cluster launch.** The substantive review was produced by backend `claude-fable-5-1`, confirmed by the client's returned model-usage metadata, in 343.33 seconds. It was one review turn with tools and MCP disabled. The CLI also recorded a small Haiku auxiliary call; the substantive review and 20,088 reported thinking tokens are attributed to Fable 5.1. This is a plan review, not independent verification of code or empirical data. The full returned text is preserved unedited below.

The earlier attempt to request broader private code/data access was rejected before submission by automatic approval review. The successful alternative shared only the self-contained utility specification, plan and overview authorized by the user. The old terminal client was too old for Fable 5.1; the already installed supported client was used in a temporary directory without granting broad home-folder trust. No project files were made available to this review process.

## Accepted and integrated

- **Utility algebra and parameter count.** The sigma=2 expression and interior-renter formula are correct. Nine initial structural coordinates plus a separately normalized psi_0 is the proposed count. The lead independently checked the expenditure constraint, first-order condition and diminishing conditional incremental material costs. Source type construction uses current dependent children when independent maturation is active (`tmp/e5f_matched_pf/code/model/intergen_eqscale_seq_optimized/solver.py:2239`).
- **The original critical path was too optimistic.** CPU-hours divided by the slot count cannot remove dependent price iterations, adaptive shock selection or final verification. Full-horizon policies cannot be promised after a late baseline. The revised plan budgets horizon/control cases explicitly inside the 12-root limit, prioritizes one candidate, allows bounded shorter screening, and separates it from certified policy output.
- **Use parallel batches inside refinements, with honest evaluation counts.** Local finite-difference matrices and candidate trial steps can run concurrently. The revised cap is 96 broad evaluations plus six refinements of at most 48 evaluations, total 384, plus the initial panel and repeats. A central nine-coordinate matrix costs 19 solves including its baseline; a claimed eight-job trial batch is not a free new Jacobian. No fixed number of optimizer iterations is promised.
- **The initial stationary pension can potentially be computed from exogenous marginals.** Structural survival fixes stationary age weights; entry income weights and the income transition are exogenous in the maintained one-market setup. This should be proved against actual marginal distributions and exploited before adding a pension root. Source: solver entry construction around 4487–4505 and `run_e5f_transition_calibration.py:624`. Actual-account verification remains mandatory.
- **Separate arithmetic balance from anticipated-income consistency.** Recomputing b from a realized population makes an arithmetic identity; households must also have anticipated that same b. The proposed test/solver reports both. A damped forward pension update is a method to test, not an a priori convergence guarantee.
- **Resolve supply consistency explicitly.** Proposed approval choice: eta=0.63 in initial and dated economies, seed-level rebasing to preserve the inherited initial price/quantity point, then recalibration of the scale. Verify units and preserve the selected supply function through all comparisons. The inherited source really does override the initial retained elasticity with 0.63 for the dated curve (`e5f_matched_pf_initial_state.py:201–209`). Rebasings at one point do not establish invariance after recalibration.
- **No-shock demographic control is useful.** It separates momentum from the preference change under the inherited observed-age initialization. Reserve an optional diagnostic root for it. It is not a reason to silently remove that initialization, and it is not needed to pretend the fitted preference shock causally explains the decline.

## Rejected or qualified

| Fable statement/recommendation | Lead assessment and resulting rule |
|---|---|
| Below-replacement fertility makes the terminal stationary endpoint impossible. | Not established for the maintained open demographic model. It includes exogenous net migration, survival and headship frozen at terminal inputs; an old positive stationary endpoint was already solved. Such affine demographic systems can have a positive finite fixed point with below-replacement births. The corrected endpoint still must be rebuilt, but do not switch to balanced growth or per-household supply on this argument. See `run_e5f_perfect_foresight_person_demography.py:315–376`. |
| The 21% model birth-count fall is already close to the correct TFR decline because female exposure grew. | Counts versus rates is the right issue, but the claimed demographic explanation was not verified. The historical model does not hold observed head-age exposure stationary. Age composition matters as well as the total number of women. Finish the female observer; do not declare the existing miss resolved or replace a target using this assertion. |
| Unchanged occupied controls with continuous saving imply a stale/cached policy bug. | Too strong. The actual saving method searches continuous intervals with piecewise continuation; an unchanged optimum at a kink can survive changes in values and slopes. Existing replay evidence does not complete certification, but this observation alone proves neither a cache bug nor its absence. Use independent anticipated-policy checks. Kernel source around 772–829 explicitly enumerates segment endpoints and interior candidates. |
| Housing cost is zero beyond the first child. | Only the additional *minimum housing requirement* is zero. The scale still increases the discretionary bundle needed for a given effective living standard, including housing. This is precisely the distinction in the author's note; do not ask him to approve a stronger claim. |
| Demote the larger-family rooms target to a check if fit is poor. | Rejected. Fit-dependent target demotion violates the identification discipline. Keep the row. Housing spending, saving, tenure and selection can respond jointly to parameters; absence of a dedicated slope does not mean the moment cannot move. |
| Center the housing floor at target response/alpha and treat smaller values as wasted. | That ratio is a conditional renter, fixed-expenditure comparison, not the empirical event-study coefficient in the full dynamic model. Keep the mapped inherited seed and broad admissible search; a value around one room is a useful additional seed, not an imposed equation or reason to discard the rest. |
| Nest the supply scale as nearly one-to-one with mean rooms. | Potential performance option only after monotonicity, units and the empirical capped-room observer are checked. Not a proven one-to-one map. Retain the scale in the nine-coordinate main specification until an equivalent profiling step is verified. |
| Never build pension derivative columns; use forward recursion. | Forward accounting must be consistent with backward expectations. Damping may work, but no contraction is established. Try it with residual checks and retain the bounded joint root as fallback; do not ban needed derivatives by assertion. |
| Set the pension consistency gate to the housing tolerance because a tighter gate can fail forever. | No general numerical theorem supports equating tolerances on different equations. Preserve separate fiscal and housing gates with explicit scaling. If convergence fails, diagnose the residuals and accuracy before changing any tolerance. |
| Replace observed-age 2007 initialization by the stationary population for tonight. | Rejected as an unrequested demographic specification change. Keep the inherited observed-age bridge and report its conditional nature. A separate stationary-population sensitivity can be proposed later; it cannot silently become the baseline to make runtime easier. |
| Finite 17-period lifetimes establish that 28 dates are adequate for historical GE. | False as an unconditional equilibrium statement. Direct individual values have finite support, but prices depend on overlapping cohorts, population and future market clearing. This distinction was already recorded in the repository. Short paths remain screening until re-cleared historical moments/prices are stable under extension. |
| One timing parameter requires age-specific rates to be identified. | Not necessarily. Multiple independent aggregate time blocks can identify amplitude and timing if the corresponding sensitivity matrix has rank. Age-specific rates are useful extra diagnostics. Identification must be checked, not inferred from a single data-type rule. |
| Calibrating early ownership necessarily estimates a bubble as owner preference. | A cyclical-confounding concern, not an established result from this review. Keep the early-window limitation and available window sensitivities explicit; do not substitute a long-run target without an approved new contract. |
| A banded price Jacobian can be recovered with 5–7 comb perturbations. | Only with justified sparsity or a controlled approximation. Dynamic demographic and price effects need not be banded. Use this only after verification; it is not an authorized default shortcut. |

## Resulting approval proposal

Approve a gated 12-hour implementation/search window, max 96 one-CPU jobs. First certify the utility/fiscal loop and empirical observer/weight contract, then the initial search; prioritize one candidate for historical shock/horizon work. Diagnostic searches and short screens must be labeled and must never overwrite the production benchmark. Retain all housing targets, observed 2007 ages and the open demographic closure; propose eta=0.63 throughout explicitly. New full policy claims require the same selected calibration, consistent expectation/fiscal paths and horizon verification.

The deliverable promised is the strongest verified result reached within the budget and a complete list of remaining failures. A promise that all calibration and policies will be fixed by morning would exceed the evidence. Source tests and numerical smoke runs occur only after the author's approval; no new model run has occurred in this review turn.

---
# Unedited Fable response

Reviewing only the plan text in your message, with no code, data, or filesystem inspected.

**Verdict.** The economics are coherent, and the utility change is smaller than its label suggests. The proposal is not launch-ready for the transition and policy stages. Stages A and B and the first calibration wave can start as declared diagnostics once three items below are settled. The wall-clock arithmetic, not core-hours, is what breaks the current plan. At old mapping times, one 100-date root alone consumes most of the night, so stage E cannot happen on the 100-date horizon.

## Blockers versus safely parallel work

Blockers, none of which can be parallelized away:

- **Terminal endpoint existence.** Flat post-2023 preferences imply below-replacement fertility, so the population shrinks without bound over the horizon. If housing supply is in absolute units with elasticity 0.63, there is no stationary rent, and the terminal rent gap sitting outside the gate is consistent with that. You must state whether supply is per household. If it is absolute, define the terminal balanced-growth path in per-household units or declare a fixed-population terminal approximation explicitly.
- **Fertility observer units.** The birth-count comparison is a units mismatch that could make you over-correct the shock. Female exposure grew over 2007 to 2023, so counts fell less than rates. A model with stationary exposure and a 21 percent count decline is already near the period TFR decline in that window. Fit rates with female exposure. Under a stationary initial population, a count target is wrong by construction.
- **Occupied-state control anomaly.** Classify it before any transition job. If saving sits on a discrete grid with argmax selection, an unchanged occupied control after a small future-income change is normal, and replacing the assertion is right. If saving is continuous, unchanged occupied controls alongside moved values point to a stale-policy or caching path keyed on occupancy. That would corrupt every transition and is a hard blocker.
- **Empirical contract certification.** Not a blocker for diagnostics, but a blocker for the calibrated label, as you already state.

Safely parallel from hour zero: adapter work, the elasticity rebase, the B panel, the first C wave as diagnostics, timing-measurement jobs, and all offline contract construction.

## Technical checks

**Utility algebra.** With sigma equal to two, the old utility is exactly minus the equivalence scale over Q, so the new specification is the old one with a different floor. Nothing about curvature changes. The renter formula is correct. Substituting the net housing bracket into a Cobb-Douglas split gives a housing share of net spending plus the floor, which simplifies to your expression. Three consequences deserve a stated decision:

- Beyond the first child, the housing cost of another child is zero, not declining. That is stronger than the sharing argument you cite. Confirm it is intended.
- The rooms gap between larger and smaller families is now generated only by selection through spending and wealth. No parameter targets it directly, so expect a poor fit there with nothing to move. Keep it, but label it a check unless the fit lands.
- The fixed-spending first-birth response equals alpha times the floor. With a nonhousing share near three quarters, hitting the reviewed room response needs a floor near one room unless spending rises at birth. Center the C wave there and treat the low end of the proposed bounds as wasted draws. If the housing choice is on a discrete grid, the objective is a step function in the floor, and the plus-minus panel will show zero or jumpy derivatives. Use a small set of discrete candidate floor values in that case.

**Parameter count.** Twelve scored rows for nine coordinates is fine on count. The supply scale is nearly one-to-one with the rooms level, so nest it inside each candidate rather than searching over it. Discount factor versus bequest strength, and first-birth fixed cost versus first-birth dispersion, are the likely weak pairs. The B panel is a finite-difference Jacobian. Compute scaled singular values from it to establish local rank, and reuse it as the initial Jacobian for Gauss-Newton refinements instead of a derivative-free search.

**Pension closure.** With a fixed tax and exogenous age-earnings profiles, the stationary pension is closed form given the age structure, and the normalization pins that structure. So the initial-state root is rent and the preference level, not rent and the pension. On the path, the pension at any date depends on fertility decided roughly five periods earlier, so recompute it by forward recursion after each mapping and iterate Gauss-Seidel. That justifies never building pension Jacobian columns. Expect the repair to raise benefits materially, which will lower private saving and push the discount factor and bequest strength up in recalibration.

| Assumption | Replacement rate of mean earnings |
|---|---|
| Old reference, 12 working to 5 retired cohorts | 0.43 |
| Balanced, if the 84 per 100 gap is all mortality | 0.51 |

Split the proposed gate in two. The identity gate is arithmetic and can sit at machine precision. The anticipation-consistency gate depends on prices cleared only to the housing tolerance, so a stricter pension gate can fail forever. Set it at the same order as the housing gate.

**Initial state.** Reweighting to the observed 2007 age distribution makes the 2007 state nonstationary. Then the no-shock baseline is itself a demographic transition, and every historical fit and policy comparison needs that control path to separate momentum from preferences. That is a second transition you have not budgeted. For overnight, use the stationary population as the 2007 state and defer reweighting.

**Early versus historical calibration.** The two-step logic is sound. Two caveats need labels. The early housing moments sit at the 2003 to 2006 price and ownership peak, so the owner preference will be calibrated to a bubble. Report a sensitivity to a long-run ownership target. Second, full anticipation front-loads the fertility drop relative to a linear preference path, which resembles the data, but the early drop coincides with the recession, so the amplitude will absorb it. With four four-year blocks, one amplitude is identified. A timing parameter needs age-specific rates, and the age gradient of the decline should be reported as a diagnostic in any case.

**Elasticity.** Use 0.63 everywhere and rebase the supply scale so the initial equilibrium point is unchanged. Both curves pass through the initial point, so the initial state is untouched and only off-path behavior changes. Do this before stage B.

## Amended plan, decisions, and fallback

Mapping time scales roughly with horizon length. Estimated per-root time at eight mappings:

| Horizon | Mapping | Root of 8 mappings |
|---|---|---|
| 100 dates | 41 to 55 min | 5.5 to 7.5 h |
| 48 dates | 20 to 27 min | 2.7 to 3.6 h |
| 28 dates | 12 to 16 min | 1.6 to 2.1 h |

Amended schedule, critical path about ten hours:

1. **Hours 0 to 2.** Stage A as proposed, plus two timing jobs running one 28-date and one 100-date mapping on the new loop with placeholder prices, plus the elasticity rebase and the anomaly classification. Hard gate at hour two.
2. **Hours 1 to 5.** B panel and the 96-case C wave together. Then refinements, but restructured. Three sequential hours at the new per-candidate time gives only about twelve evaluations per chain, which is below one simplex in nine dimensions. Run each chain as batched Gauss-Newton with eight concurrent trial jobs per iteration, so twelve chains use the same 96 slots and finish four to five iterations.
3. **Hours 5 to 8.** Stage D at 28 dates only, twelve roots in parallel, pension by recursion, cap eight mappings. Since maximum lifetime is about seventeen periods, households deciding through 2023 do not care past the 2090s, so 28 dates is adequate for historical fit and the horizon check is the confirmation.
4. **Hours 7 to 11.** Pick the best one or two. Launch their 48- and 100-date roots warm-started from the 28-date path as the horizon check. In parallel, launch the property-tax baseline and reform at 28 dates, labelled screening, valid for roughly 2023 to 2060 given terminal contamination. Supply and LTV variants only if slots and time remain.
5. **Hour 11 to 12.** Verification and readout.

Checkpoint every mapping so a preempted job restarts mid-root. The queue is empty but nothing is reserved, so treat preemption as expected. For the price Jacobian, estimate a banded matrix with a comb of simultaneous perturbations in five to seven mappings, then Broyden updates.

Decisions for you to settle before launch:

- Stationary 2007 population for overnight, with reweighting deferred.
- Terminal endpoint units, and whether supply is per household.
- Rate-based fertility targets with female exposure, counts as a check.
- 28-date screening for D and E overnight, with long horizons only for finalists, and policy results labelled screening.
- Elasticity 0.63 everywhere with a rebased supply scale.
- The rooms-gap moment as a check unless the fit lands.
- Zero marginal housing cost beyond the first child as an intended feature.
- The anticipation test replacement, conditional on the grid classification.
- Acceptance of bubble-peak ownership as the target, with a stated sensitivity.

Fallback if the contract is not certified: run B, the C wave, and the refinements against the currently available moments with declared absolute tolerances, labelled diagnostics. Run D as a diagnostic amplitude fit to published period TFR by four-year block with a declared tolerance. Do not launch E. The morning product is then the certified fiscal loop, the timing table, the local rank diagnostic from the B panel, and the horizon comparison. Those are prerequisites regardless of the contract, so nothing is wasted.

Bottom line: the utility and pension repairs are coherent, the 100-date transition budget is infeasible overnight, and three definitional questions, the terminal endpoint, the fertility units, and the occupied-state anomaly, must be answered before any transition job starts.
