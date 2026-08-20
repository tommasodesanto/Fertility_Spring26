# Active Theory and Quantitative Decision Ledger

Updated: 2026-08-19

This is the short working ledger for the current paper. It records decisions
that can still change the theory, quantitative interpretation, or paper-facing
claims. `CALIBRATION_STATUS.md` remains the source for numerical results and
reproducibility contracts.

## Active points

1. **Child unit resolved.** The live sequential model's parity states count
   literal children and match CPS children-ever-born per woman. The transition
   therefore converts top-code-adjusted child births into future entrant
   households by dividing by `2.1`. In the toy economy, one household-child unit
   corresponds to `2.1` quantitative-model children. The divisor is an external
   replacement normalization, not an explicit childhood-survival state.

2. **Finish the small transition theory.** Write the explicit population law,
   the impact response to lower desired fertility, demographic momentum, the
   housing-price feedback, and the condition under which no positive closed
   steady state exists. The diagrams must follow that law rather than remain
   schematic.

3. **Validate current fertility flows.** Construct the exact model counterpart
   to period TFR and assess whether the model makes post-2023 fertility too low.
   Completed fertility for an older cohort cannot by itself identify the speed
   of the future population transition.

4. **Address the two binding housing-quantity misses.** The first-birth housing
   response and mean rooms remain the largest contributions to the current
   calibration loss. A fresh 23-task radius-`0.01` coordinate panel found no
   material local improvement: the best one-coordinate loss is `36.0222`
   versus `36.0992`, and the eleven-column Jacobian has rank ten at the
   predeclared `1e-3` threshold. Its weakest direction is primarily the already
   low bequest shift. Do not force a ridge update through this weak direction.
   The next structural test should separate the first-child room increment from
   the level of large-home demand or from the tenure utility margin; any added
   parameter requires an identifying moment and a fresh target/code contract.

5. **Quantitative presentation.** The current 2023-forward baseline
   uses household totals and age shares as imposed dated inputs through 2023.
   The fixed four-group alternative is rejected as a demographic closure: it
   resets all age-cell masses, removes most raw youngest entrants, and identifies
   household-formation/headship wedges rather than immigration. The additive
   semi-open prototype preserves endogenous entrants but reveals that an
   explicit household-formation margin is needed. The main exercise is the
   fitted-2023 closed finite-horizon path; the open path is a sensitivity.

6. **Choose and state the long-run interpretation.** Neither the fitted-2023
   baseline nor the four-group historical calibration has a verified positive
   closed terminal steady state: their audited maximum renewal ratios are
   respectively 0.7085 and 0.7094. The closed national paths are finite-horizon
   contraction benchmarks. The open paths are sensitivities to an externally
   normalized entrant flow. Do not describe either as something stronger.

7. **Update the paper and slides after Points 1--3.** Present the old steady
   state, impact response, demographic momentum, and endogenous housing
   feedback in that order. Do not claim movement between two positive steady
   states unless a positive closed root is actually established.

8. **Keep structural policy analysis separate.** The rebated property-tax
   transition may be shown from 2023 through 2035. Lead with the maintained
   national supply elasticity `1.75`; the Coven California elasticity `0.232`
   and the Baum-Snow--Han U.S. urban elasticity `0.5` are impact sensitivities.
   The path is a finite transition with equal rebates, not a stationary or
   welfare comparison.

## Legacy points -- to be reviewed

These were previously recorded as unresolved. They are not automatically part
of the current specification; each must be either revalidated against the live
model or retired explicitly.

1. **Outside-entry anchor.** Decide whether the outside-origin object is an
   age-aggregated net-migration equivalent or a literal age-18 entrant flow.
   The old `0.169` value is presently only an old-state normalization, not a
   constant realized migration share.

2. **Open-closure counterfactual contract.** Reassess the older assumptions of
   policy-invariant outside inflow and policy-invariant local-born retention
   before using an open closure for policy.

3. **Balanced-growth closure.** The older proposal to use a balanced-growth
   path was left unresolved because level housing supply, proportional supply,
   and Stone--Geary housing floors do not fit together innocuously. Review only
   if the current contraction interpretation is rejected.

4. **Housing-supply normalization.** Reconcile the one-intercept representation
   `H^S = H_0 p^xi` with older paper, slide, and packet presentations that use a
   separate reference price. This is a reporting convention, not a reason to
   re-solve the model.

5. **Entry-response primitives.** The older entry taste scale and local-born
   retention weight were not empirically disciplined. They must not return as
   production parameters without targets or explicit external restrictions.

6. **Geographic interpretation.** Separate national demographic renewal from
   local spatial reallocation. A city entrant flow or migration response cannot
   be relabeled as a national population mechanism.

7. **Funded policy and purchase-grant architecture.** The previous funded
   workflow remains disabled. Its renewal accounting must use the live unit
   convention, and purchase-grant eligibility requires a one-time/past-owner
   state to prevent cycling.

8. **Preserved one-shot model.** Keep the circulated one-shot model as a
   documented fallback and comparison, not as the source of current sequential
   transition results.

9. **Older wealth, bequest, and numerical-boundary concerns.** Reopen only if
   they affect the live calibration diagnostics. Historical loss values and
   conclusions from superseded target systems are not comparable to the current
   transition calibration.
