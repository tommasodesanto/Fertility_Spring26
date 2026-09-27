# Same-parameter borrowing diagnostic

Reference: authenticated overnight continuation de_0093. All economic
parameters, including psi_child=0.14416490417555738, remain fixed. The original
14-row target system remains a reporting comparison, not an objective being
minimized here. Fertility is not renormalized. No change to entry wealth,
preferences, prices imposed externally, housing grids, target weights or supply.
Equilibrium prices respond endogenously.

The experimental arm removes both purchase screens and the collateral/unsecured
saving floors. It replaces them with continuation feasibility for every
positive-probability income/child outcome and nonnegative liquidation wealth
at every possible death: b_next + (1-selling_cost)*price*house >= 0.
Bequest preferences retain their existing gross valuation; the net valuation
here concerns debt repayment and the estate funding ledger. No default or
insurance mechanism is introduced. This repayment requirement can itself
restrict unsecured debt, so it is not an unconditional relaxation at all states.

The boundary is the first feasible wealth-grid node. The native value cutoff
V > -1e9 classifies feasibility, potentially rejecting very negative finite
utility. This is a conservative discrete diagnostic, not a certified continuum
frictionless optimum. Impossible intervals halt rather than being smoothed over.
Birth/entry imbalance and estate funding shortfalls are reported explicitly;
no closed-demographic-renewal counterfactual is claimed.

Three serial cases: reference replay, experimental equilibrium, exact repeat.
Full operator, household-budget, probabilities, market and pension checks remain.
The collateral audit is replaced by an independent exact transaction and estate
solvency audit. Full14-target and31-parameter tables and17standard figures are
required for successful cases. Eleven focused adapter tests pass.

`run_v1` preserved: reference equilibrium solved in56.41seconds, market residual
2.12e-8; all31 parameter estimates exactly equal the reference and all14 model
moments within2.34e-6. Wrapper rejected the runtime eq_iter counter as a changed
parameter after completing model checks. This is a workflow error, not a failed
equilibrium. Its original runner is retained as runner_snapshot.py.

`run_v2` excludes only eq_iter from the all-parameter-field assertion; solver.py
sets it on each price evaluation. No economic change. Original absolute deadline
1790522449.449255 retained (11:20:49 EDT), one worker,10minutepercase. plan.json
pins the reviewed runner/helper and reference. Check checkpoint.json and
complete.json for actual completion; do not infer results from launch.
