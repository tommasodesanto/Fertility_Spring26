# Initial stationary PAYGO review — September 11, 2026

**Verdict: source and pure-check PASS for the maintained one-market initial stationary solver.** The age/income shortcut and supply rebase are correct. This review does not certify a solved equilibrium: the lead owns the saved-distribution marginal audit and exact-candidate cluster smoke. No model solve, source edit or cluster launch was performed by this reviewer.

The main status and approved plan were read, including the September 11 approval that supersedes the older preparation hold. Reviewed code is the isolated `tmp/e5f_matched_pf` checkout. The observed-age announcement state and terminal person-demographic closure remain outside this shortcut's authorized use.

Reviewed final SHA-256: `e5f_stationary_paygo.py` = `d723bae3e2a6d912f1095a5aa94a1780eb9ad6a43babdd45659bea31c680a6f9`; `test_e5f_stationary_paygo.py` = `7d7508c85991b07e4351ac95c22a0b2444447afa7ce8971adb392e23ac237f5f`.

## Why the marginal recursion is valid

Let \(a_j\) be the row vector of household mass at age \(j\) over earnings states, before a common normalization. The actual forward code implies

\[
a_0=E\pi_0,\qquad a_{j+1}=s_j a_j\Pi_z,
\]

where \(E\) is one-market entrant mass, \(\pi_0\) is normalized `z_weights`, and \(s_j\) is the scalar age-survival probability. The shortcut sets \(E=1\), propagates this recursion, and divides by total mass. This gives the same normalized marginal.

Source evidence in `intergen_eqscale_seq_optimized/solver.py`:

- Entry uses `entry_by_loc[i] * z_weights[zz]` with wealth weights summing to one (4494–4509); income-dependent entry wealth changes wealth allocation, not total entry mass by income.
- `_censor_entry_dead_mass` (3947–3974) relocates only the wealth index within each fixed tenure/location/income/parity/child column. It never drops or changes income mass; columns with no feasible frontier remain subject to the existing rejection gate.
- Fertility redistributes a parent's mass among parity/child states. The forward survival factor is age-only (4891); location, tenure and saving moves preserve mass. With one location there is no spatial selection of the tax base.
- Stochastic child maturation preserves parent mass; the independent-count transition is binomial and row-stochastic, including self-loops for unreachable cells (`parameters.py`, 896–928). Mature children are separately counted, not injected into this stationary head recursion.
- `Pi_z` is applied on the age transition (4980–5023). The recursion correctly handles entrants whose income weights are not invariant under `Pi_z`; its test uses alternating earnings states. Active parameter construction ordinarily enforces invariant weights, but that stronger assumption is unnecessary here.
- Aggregate normalization (5053–5069) is a common scalar. Current housing-choice realization (5071–5088; 3619–3653) preserves age/income cells in one market. Therefore both beginning and realized current distributions have the same relevant marginal. Tiny mass-pruning thresholds are numerical exceptions, which the actual-marginal gate appropriately checks rather than assuming exact floating-point equality.

## Pension units and binding

`fiscal_accounts` integrates worker earnings in period units, \(\Delta w a_j z\), and retiree exposure \(1+s_z(z-1)\). It therefore returns the period benefit

\[
b=\frac{\tau\Delta\sum_{j<J_R,z}w a_j^{\rm income}z\,m_{jz}}
{\sum_{j\ge J_R,z}[1+s_z(z-1)]m_{jz}}.
\]

Common entrant/population scale cancels. The binder sets working income to \(\Delta(1-\tau)wa_j^{\rm income}\) and retired base income directly to \(b\). `income_at_state` (228–235) then applies the correct earnings/retirement multiplier. There is no second multiplication of pension by four.

The direct `solve_markov_income_equilibrium` route (1164–1283), its shared precomputation and full household solve do not rerun the annual pension resolver or `apply_overrides`; the bound income survives to actual household optimization. Calling the general parameter constructor after binding would be unsafe. The new post-solve income-array, location-pension and period-marker check now rejects that inconsistency explicitly.

## Supply units and findings resolved during review

The actual fast/full solution packers use \(H=H_0[(q+\delta+\tau_H)p/\bar r]^{\xi}\) (7143–7144 and 7177–7178). Thus the implemented mapping
\(H_0^{new}=H_0^{old}R^{\xi_{old}-.63}\), with \(R=\texttt{user_cost_rate}\,p/\bar r\), preserves the old stock at the supplied asset price. It preserves one point, not the whole supply curve. The test was independently compared to the actual fast solution packer.

Two guard gaps were reported and fixed by the lead during review: the certificate originally did not check the actual income array, and the wrapper accepted any supplied market tolerance. The updated code verifies binding consistency and requires `0 < tol_eq <= 2.5e-5`. The corrupted-income regression now fails as intended.

## Verification and remaining boundary

All four revised `test_e5f_stationary_paygo` tests pass with JIT disabled. Independent pure helper checks also passed: common-mass rescaling by 0.03, 1 and 1,000,000; rejection after multiplying retired household income by four; exact preservation of income-cell entry mass by the actual censor helper; and rebased supply equality against the actual solver packer. No household optimization ran.

The helper test suite does not yet exercise `solve_balanced_initial_equilibrium` through a fake solver (correct route, false convergence and loose-tolerance rejection). That would be useful regression coverage; the more important pending evidence is the lead's exact-candidate smoke, with actual marginal, fiscal, market, feasibility and existing household gates. The shortcut must not be reused after empirical age reweighting or in a dated population path without recomputing the actual fiscal account.
