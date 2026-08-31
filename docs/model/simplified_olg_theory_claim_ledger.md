# Simplified OLG theory claim ledger

This ledger records exactly what the paper-facing two-generation theory proves,
what is conditional, and what is established only by the reproducible mixed-
tenure construction. The source files are
`latex/simplified_olg_paper_theory_section.tex` and
`latex/simplified_olg_paper_theory_appendix.tex`; the numerical driver is
`code/model/tools/build_simplified_olg_mixed_tenure_theory.py`. The independent
audit and its resolution are recorded in
`docs/model/simplified_olg_theory_hostile_audit_20260831.md`.

| Claim | Status | Conditions | Paper location | Verification |
|---|---|---|---|---|
| The old renter and old owner problems have the stated active-set solutions. | Proved | Positive user costs and resources; the row-specific complementary-slackness inequalities. | Lemma `lem:olg_old_allocations` and equations `olg_app_old_*` | Direct substitution into the KKT system. |
| Interior saving collapses each conditional young problem to the lifetime budget `(1+k)c + chi n + r_m h = W`. | Proved | Young saving and the relevant next-period old constraints are interior. | Proposition `prop:olg_app_conditional_policies` | Mortgage inflow and repayment cancel exactly; saving policies are reported separately for renters and owners. |
| Unconstrained conditional consumption, fertility, and housing have the closed forms in the paper. | Proved | Interior saving, old choices, and tenure-specific housing limit. | Equation `eq:olg_unconstrained_policies` | Algebraic solution of the three first-order conditions. |
| Fertility at a binding housing limit is unique. | Proved | Nonempty log domain and positive user cost. | Equation `eq:olg_binding_fertility`; proof of Proposition 3 | The scalar first-order condition falls strictly from positive to negative infinity over its feasible interval. |
| The logit tenure formulation produces positive renter and owner shares and continuous aggregate choices. | Proved | Finite conditional values and taste scale `kappa_T>0`. | Equation `eq:olg_owner_share`; Appendix A.2 | The logit probability is strictly inside `(0,1)`; conditional policies are continuous within their stated active sets. |
| Every positive steady state satisfies average fertility `1/nu`, has zero taxable real gains, and solves the two-equation `(P,T)` system. | Proved | Positive stationary cohort masses and stationary household policies/distribution. | Proposition `prop:olg_steady_state` | Follows from the two cohort laws, housing clearing, and the balanced government budget. |
| A positive steady state exists under the stated fiscal contraction and replacement-boundary inequalities. | Conditional theorem | Compact feasible price-transfer rectangle, continuous aggregate policies, fiscal self-map/contraction, and a replacement sign change. | Proposition `prop:olg_app_ss_existence` | Banach fixed point for the fiscal locus plus the intermediate value theorem. These conditions are sufficient, not automatic. |
| A positive steady state is locally unique and differentiable in the fertility preference when its two-by-two Jacobian is nonsingular. | Conditional theorem | Differentiability at the root and `det J != 0`. | Equation `eq:olg_app_ss_jacobian` and comparative static `eq:olg_app_ss_comparative_static` | Implicit-function theorem. No unconditional sign is claimed for the price or population response. |
| A terminally closed finite-horizon truncation exists locally for a small one-sided permanent shock on a selected appreciation branch. | Conditional theorem | Regular endpoint, strict household active sets, nonsingular branch Jacobian, and strict directional price signs consistent with every selected gains branch. | Proposition `prop:olg_app_local_transition` | Directional branch construction plus the implicit-function theorem. This is a root through the terminal date, not an exact infinite-horizon transition. |
| A terminally closed finite-horizon truncation exists globally under coordinate-wise boundary signs. | Conditional theorem | Continuous stacked residual on a compact feasible box and opposite signs on each pair of coordinate faces. | Proposition `prop:olg_app_global_transition` | Poincare--Miranda theorem. The boundary signs must be verified for a parameterization. |
| A pointwise limit of terminally closed truncations is an infinite-horizon transition equilibrium under compactness, tightness, fixed-date convergence, continuity, and uniform convergence of the finite-horizon tails to the terminal steady state. | Conditional limiting theorem | The convergence conditions and uniform-tail equation stated in Appendix A.4. | Proposition `prop:olg_app_transition_limit` | Fixed-date passage to the limit plus the uniform-tail condition. A terminal condition at a moving horizon alone is not sufficient. |
| One permanent shock can generate taxable gains for several cohorts. | Proved, conditional on the equilibrium price path | `P_s>P_{s-1}` on each stated date. | Corollary `cor:olg_app_one_shock_gains` | Immediate from `ell_s=tau_g(P_s-P_{s-1})`. A two-period life determines the basis, not the duration of price adjustment. |
| The current-goods marginal-value gap between a constrained young owner and an interior old incumbent is `zeta_OF + ell_t + q ell_{t+1}`. | Proved | Binding young down-payment constraint, slack owner-size limit, interior young saving/future old choices, and interior incumbent retention. | Equation `eq:olg_marginal_values` | Young and old housing KKT conditions. |
| A small intergenerational housing transfer is a compensated improvement when the marginal-value gap exceeds the real transfer cost. | Conditional welfare proposition | Fixed fertility, cohort masses, tenure, and unaffected allocations; pledgeable bridge finance; two-date lump-sum transfers; incremental rebate recovery; marginal resale next period. | Proposition `prop:olg_reallocation` and Assumption `ass:olg_two_date_ledger` | Willingness-to-pay/willingness-to-accept argument with an explicit two-date government ledger. |
| Relaxing a binding housing limit raises fertility exactly under the reported inequality. | Proved partial-equilibrium result | Tenure, prices, lifetime wealth, and active set fixed; interior fertility. | Equations `eq:olg_fertility_derivative` and `eq:olg_fertility_sign` | Implicit differentiation; analytic derivative agrees with a centered finite difference to `8.0e-12` in the construction. |
| The displayed renter-owner path is a high-accuracy terminally closed truncation of the specialized numerical model. | Constructively verified | Parameter table in Appendix A.7; final steady-state prices, transfers, and continuation values after date 28, with the induced terminal state recorded rather than imposed. | Figure `fig:olg_mixed_transition` and Appendix A.7 | Max scaled residual `2.07e-10`; terminal state gap `1.10e-8`; finite-root Jacobian minimum singular value `0.401`; 24-vs-28 first-nine-date gap `1.84e-15`. These diagnostics do not verify the uniform-tail theorem. |
| Both tenures and all asserted active sets remain strict in the numerical truncation. | Constructively verified | Same parameterization. | Appendix A.7 | Owner share `0.4704`--`0.7877`; positive saving; strictly binding young caps; owner-size slack at least `1.2046`; positive old-renter cap wedge; positive old-owner sales, user cost, and estate slack. |

## Claims deliberately not made

- No global welfare ranking is made for policies that change fertility or the
  number of future households. Such a ranking requires an explicit
  endogenous-population social criterion.
- No claim is made that a positive steady state or a transition exists for
  every parameterization, that either object is globally unique, or that house
  prices must rise after every fertility shock.
- The computed horizon-28 path is not labeled an exact infinite-horizon
  equilibrium. It is a terminally closed approximation with small residuals,
  a small terminal-state gap, and horizon stability over the reported window.
- The partial-equilibrium fertility derivative does not sign a funded policy
  in general equilibrium, because prices, transfers, tenure, and lifetime
  resources can also move.
- The mixed-tenure construction is a theoretical existence and mechanism
  example, not an empirical calibration.
- Nothing in this simplified construction establishes a positive closed
  stationary root for the active quantitative lifecycle model. Its live
  stationary and transition status remains governed by `CALIBRATION_STATUS.md`.
