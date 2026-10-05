The v5 failure is a date-0 inherited-state rejection, not an OOM or root-budget failure.

1. Failing operation and mechanism

- `map_001` fails in `_require_exact_inherited_distribution`, which rejects occupied cells when \(V\le -10^9\) and their mass exceeds \(10^{-12}\): [run_dynamic_population_transition.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/run_dynamic_population_transition.py:415)–[453](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/run_dynamic_population_transition.py:453).
- The receipt records 2,823 affected cells and \(9.96425726572\times10^{-5}\) mass—about \(10^8\) times the gate—not projected or altered: [inherited evidence](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/transition_readiness_v1/current_baseline_20261003/two_shock_v1/execution_smoke_v5/results/fit_fit_v1/run/stage1/candidate_0001/horizon_024/map_001/inherited_state_evidence/inherited_infeasible_1791157483239129132.json:1). `failure.json` confirms this is terminal under the exact-mass rule.
- It is date 0 (2007), not a later propagated distribution: the inherited state is copied before the forward loop, which begins at `period=0`; the policy is then evaluated against that untouched `state.g_pre` before any forward advance: [perfect-foresight runtime](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/transition_readiness_v1/current_baseline_20261003/two_shock_v1/execution_smoke_v5/frozen/source/code/model/tools/run_e5f_perfect_foresight_transition.py:675)–[713](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/transition_readiness_v1/current_baseline_20261003/two_shock_v1/execution_smoke_v5/frozen/source/code/model/tools/run_e5f_perfect_foresight_transition.py:713).
- V5 uses the stationary endpoint constants for every cold-path date: [two_shock_runtime.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/transition_readiness_v1/current_baseline_20261003/two_shock_v1/execution_smoke_v5/frozen/source/code/model/experiments/birth_count_choice/two_shock_runtime.py:101)–[107](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/transition_readiness_v1/current_baseline_20261003/two_shock_v1/execution_smoke_v5/frozen/source/code/model/experiments/birth_count_choice/two_shock_runtime.py:107): \(q_t=0.6482672107\), pension \(=0.9177840475\). This changes the date-0 Bellman continuation problem, while preserving the inherited population exactly.

2. What the negative values establish

They do not establish a current-budget violation. The gate only observes \(V\le-10^9\), not consumption or the dated budget audit. The Bellman continuation explicitly makes a state infeasible if any reachable continuation or death branch is infeasible: [household.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/transition_readiness_v1/current_baseline_20261003/two_shock_v1/execution_smoke_v5/frozen/source/code/model/experiments/birth_count_choice/model/engine/household.py:654)–[674](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/transition_readiness_v1/current_baseline_20261003/two_shock_v1/execution_smoke_v5/frozen/source/code/model/experiments/birth_count_choice/model/engine/household.py:654).

Thus the large finite negatives are consistent with propagated future infeasibility. No date-0 budget audit was written because rejection occurs before fertility/current-choice/audit execution. The census’ indebted-owner concentration is evidence that the endpoint-flat forecast breaks their dated continuation support; it is not proof of a contemporaneous owner-budget violation.

3. Smallest numerical continuation recommendation

No proven feasible initializer can be derived from these artifacts alone. In particular, the renter condition depends on adjacent prices:
\[
r_t=u q_t+(q_t-q_{t+1}),
\]
not only on \(q_t\): [perfect-foresight runtime](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/transition_readiness_v1/current_baseline_20261003/two_shock_v1/execution_smoke_v5/frozen/source/code/model/tools/run_e5f_perfect_foresight_transition.py:358)–[386](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/transition_readiness_v1/current_baseline_20261003/two_shock_v1/execution_smoke_v5/frozen/source/code/model/tools/run_e5f_perfect_foresight_transition.py:358). The v4 affine failure already shows that merely restoring \(q_0=0.7761012531\) is insufficient for renters.

The smallest non-economic-change hypothesis is therefore a date-0 continuity initializer, not another endpoint-flat path: retain the reference initial price and pension, set \(q_1=q_0\) so date-0 rent equals its stationary value, then smoothly bridge \((q_t,b_t)\) to the solved endpoint. This is only a feasible-initializer hypothesis: it preserves the renter’s date-0 rent identity, but no available arithmetic proves inherited-owner continuation feasibility or later-date renter feasibility.

4. Exact check before any empirical retry

Run exactly one unrooted, 24-date native mapping at the initial stage-\(1\) \(\psi\), with the immutable sources/gates and exact inherited state, using that continuity path. Do not run the scalar fitter or root iterations. Require:

- zero inherited infeasible mass at the existing \(10^{-12}\) threshold;
- all dated household-budget, purchase, estate, and projection gates;
- saved first-failure evidence if it rejects.

This one mapping is the minimum native validation needed to distinguish a viable numerical initializer from an unsupported hypothesis. If it fails, its date-specific evidence—not a mass cutoff, projection, or gate change—must determine the next continuation adjustment.

Unknown: the exact inherited-owner feasibility boundary in \((q,b)\)-path space; later-date renter budget margins; and whether any continuous bridge reaches the endpoint within 24 dates while retaining every current gate.