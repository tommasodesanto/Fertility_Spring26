# Terminal steady states and a clean shock experiment

September 13, 2026. Read-only methods assessment; no new numerical job, calibration change, or package migration was launched for this review. The author's 12:30 EDT technical cutoff remains in force.

The proposed benchmark is a transition from a verified initial steady state after one unexpected, permanent change in the fertility preference. This removes estimation of the four historical shocks from the experiment. A temporary shock can instead return to the original steady state; a permanent change generally requires a different terminal steady state.

## Boundary and algorithm

First solve the terminal equilibrium at the declared permanent preference, with the same household, demographic, housing-supply, PAYGO and property-tax-rebate rules as the transition. Set its lifetime value functions and terminal aggregate conditions as the far-future boundary. Guess all intervening prices, pensions and rebates; solve households backward and carry the inherited population forward; update the paths until dated markets and budgets clear. Inspect the actual carried endpoint against the terminal equilibrium and compare the early response across horizons. Do not overwrite the carried distribution with the terminal distribution.

This is a boundary-value problem solved by time-path iteration or a root method over the whole path. Classical shooting instead adjusts initial jump variables and integrates forward. Both use boundary information, but their numerical implementations differ.

The current driver, `code/model/tools/run_e5f_final_rebated_history.py`, instead solves the endpoint price, pension and rebate jointly with the path and constructs lifetime policies at those constant conditions. It clears the actual endpoint but does not require the population to reproduce itself afterward. Increasing its horizon to 100 does not, by itself, replace that boundary with a terminal steady state.

## Findings from the live implementation

- The historical initial-state builder reweights the stationary household distribution to observed 2007 ages. A clean stationary shock experiment must omit historical reweighting and conditioning and verify a no-shock path under the same population law used in the experiment. The existence of a stationary household calibration does not establish that the full demographic state is stationary under that law.
- The retained `e5f_balanced_terminal.py` adapter inspected in the frozen `corrected_history_source_v2` explicitly requires zero property-tax rebate. The live baseline rebates the tax. The canonical status also records that its older fixed-migration terminal construction cannot be reused unchanged for the zero-migration baseline. These are implementation and economic-closure issues, not reasons to relax numerical tolerances.
- With zero migration, a positive stationary population at the new preference must be established. It is not valid to restore migration, change the preference to force replacement, or normalize away population decline merely to obtain an endpoint. Failure to find a root is also not proof that no root exists.
- Saved native receipts on September 13 show 100-period full mappings taking 1,735.79 and 1,972.77 seconds, with 200 dated household solves. A completed 24-period fixed-preference root took 5,215.16 seconds over 16 evaluations. These are observations of the current implementation, not lower bounds for a better algorithm. A verified new 100-period transition cannot responsibly be promised before 12:30 today.

## Relevant primary sources

[Auclert, Bardóczy, Rognlie and Straub (2021)](https://straub.scholars.harvard.edu/file_url/106) develop sequence-space Jacobians and use them for linear impulse responses and nonlinear perfect-foresight transitions. Section 6, equation (38), updates the entire path using a steady-state Jacobian. The subsection on a transition to a new steady state uses the terminal-state Jacobian for a permanent shock. Its speed requires constructing the appropriate model derivatives.

The official [sequence-jacobian package](https://github.com/shade-econ/sequence-jacobian) provides both linear and nonlinear routines. In [Block.solve_impulse_nonlinear](https://github.com/shade-econ/sequence-jacobian/blob/master/src/sequence_jacobian/blocks/block.py), lines 149–187 as inspected today, the Jacobian is computed or supplied at `ss`, factored, and reused while nonlinear residuals are reevaluated. The separate `ss_initial` argument supports different initial conditions. This solver does not establish the economic validity of an externally supplied terminal steady state.

The independent subagent checked the paper's temporary and permanent shock distinction against the package's lead/lag padding and heterogeneous-household initialization. The lead independently inspected the nonlinear whole-path update in the official source and checked our live terminal and initial-state adapters. No package benchmark was run. Reusing the terminal-state Jacobian to update our nonlinear path is a plausible development direction, not a verified speedup; affordability thresholds and discrete choices also need derivative and convergence checks.

Our model would still need a verified representation of its household decisions, endogenous births, population masses, entry rules and fiscal accounts. A package installation alone cannot supply those objects. The authors' [life-cycle extension code](https://github.com/Mv77/LC-SSJ_public) is a relevant follow-up source; it has not been integrated or benchmarked against this model.

## Recommended first experiment

Verify a no-shock path from a fully consistent initial steady state, then solve one declared preference shock with an equally consistent terminal equilibrium. Preserve the original parameters, supply elasticity, taxes and numerical gates. Report fertility, population, household counts, housing services, prices, rents, pensions and rebates together with equilibrium residuals and endpoint distances. Attempt a short and a 100-period version only after the endpoint and exact-loop checks pass. This is a mechanism and solver test, not a historical shock fit or a replacement for it.
