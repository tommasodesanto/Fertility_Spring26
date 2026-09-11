# Anticipation test repair — reviewed implementation, cluster verification pending

Owned source: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_matched_pf/code/model/tools/test_e5f_social_security_compiled.py`.
SHA-256: `414f7d3e1320febd69809fd4114a94d649dbdb7794b24c1127b4eebc6a8908ab`.

The recovered CLI-worker draft was independently reviewed and strengthened. It had verified the outer continuation handoff, but its saving-interiority test checked only the wealth-grid endpoints and provided no conditional tax-control acceptance check. The revised suite has seven test methods sharing the unchanged six two-date paths: 24 dated Bellman calls plus the original tiny-fixture initialization. Observation wrappers call the actual compiled kernels; they do not substitute policy values or add solves.

## Precisely named acceptance tests

`test_dated_fiscal_income_and_age_support_reach_both_bellman_passes` checks:

- Four actual calls in backward dates 1/0 and forward dates 0/1, including exact value replay and the common terminal value.
- Independently reconstructed period pension/payroll income and its household income-state multiplier inside the actual compiled resource vectors. This observes the two working probes plus a declared age-5 pension recipient in every call.
- Identical current income across each shock/control pair and exactly dated continuation handoffs.
- The one-step age support implied by the two-date experiment: a future pension affects current age index 4, and a future payroll tax affects indices 0–3. Unaffected ages remain unchanged; pension values weakly increase and payroll-tax values weakly decrease, with strict responses somewhere and at occupied current states.

`test_predeclared_controls_satisfy_exact_continuation_and_crossed_policy_incentives` checks:

- Pension state `(1,0,0,4,0,1,0)`, previously saving `0.9382044371 -> 0.5013730264`; tax state `(10,1,0,0,1,0,0)`, previously saving `0 -> 0.8407770257`. Axes are wealth index, conditional tenure, location, age index, income state, parity, dependent-child count. These states were chosen using the previous diagnostic, before this new run; they are not selected by a maximum in the new output.
- Actual kernel continuation equals the independent expectation
  \(\bar V(b')=\sum_{z',m'}\Pi_z(z,z')\Pi_m(m,m';n)V_{t+1}(b',\tau,i,j+1,z',n,m')\)
  in both current-date Bellman passes. The tiny fixture has no survival mixture or readiness gate; that fixture contract is explicit and fails if changed.
- Recorded conditional values equal independently reconstructed \(u(c,h)+\beta\bar V(b')\), reported consumption/housing match the objective bundles, and report-only consumption/housing floors are inactive at these probes.
- Actual saving bounds come from current resources, housing/nonhousing floors, collateral and debt constraints, not just grid endpoints. Pension saving is interior under both paths. Tax saving releases a verified lower-bound corner into the interior; the test does not claim it was initially interior.
- All current kernel arguments and feasible-set inputs are identical across each pair except the continuation array. Each selected policy strictly beats the other policy under its own continuation. At fixed saving, the entire value difference is exactly \(\beta\Delta\bar V\); at the changed saving, the full difference decomposes into current utility and discounted continuation changes.

Existing explicit-baseline neutrality, first-date fiscal balance, household budget, population mass, replay, and entry-queue gates remain intact. No solver, utility, target, weight, population law or fiscal equation was edited.

## Preserved failures and verification

Original smoke `17352615` remains failed: its occupied-control maxima were zero for pensions and roundoff for taxes. Its source is retained at `a654219c`. Original error SHA-256 is `7d8384c1b857e7539712e2848a1f390e74c2c6e991974125fc79b60e1161e181`. Diagnostic `17353361` remains unchanged; JSON SHA-256 is `22811cadf2a906d1df39e7767c6ee4e9a21e8be8bcf5d23951a00528910f0420`. The source docstring points explicitly to both historical jobs and their preserved main-checkout evidence directory.

Passed locally: final-source syntax compilation and scoped whitespace check; ten AST-isolated scalar checks of renter/owner utility, capped/uncapped renter housing, budget bounds, constant-continuation shifts, and invalid-state/infeasible-consumption rejection. The helper was extracted without importing the model or constructing the fixture. Initial system-Python bytecode writing encountered a cache-directory permission error; compilation to an explicit temporary output and the final in-memory compilation both passed.

**Not yet certified:** no local compiled/model run or cluster launch occurred. Lead must inspect the diff and run this exact six-case smoke on Torch. New crossed-policy gain and monotone-value checks require real numerical verification; if one fails, preserve and diagnose it rather than lowering thresholds. This is the unchanged prior tiny Stone–Geary fixture, so it checks fiscal anticipation and accounting, not the separately revised parenthood equivalence-scale utility. It proves neither an occupied saving response nor general global optimizer optimality, equilibrium, or policy validity.
