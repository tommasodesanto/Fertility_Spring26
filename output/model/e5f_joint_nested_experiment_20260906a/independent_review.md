Stopped on an unclosed implementation assumption.

1. **Location choice is omitted from the four-plan definition.** `q(\tau,x)` is only a housing-product maximum ([experiment lines 16–21](file:///Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/docs/model/e5f_joint_nested_experiment.md#L16)), but the live solver constructs tenure values and then applies a destination-location logit ([solver lines 2622–2689](file:///Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/intergen_eqscale_seq_optimized/solver.py#L2622), [2691–2717](file:///Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/intergen_eqscale_seq_optimized/solver.py#L2691)). Before implementation, specify whether \(q\) includes a tenure-conditional location inclusive value (preserving the existing location shock after the four-plan choice), whether location is jointly chosen, or whether it is fixed. Otherwise it is unclear how four alternatives reproduce the intended menu.

2. **The “reversed sequential-revelation” control needs an explicit formula and label.** Lines 85–89 say it uses the same \(Q\)’s with tenure scale \(\kappa\) and fertility scale \(\lambda\kappa\), but do not define it. If intended as fertility-first sequential logit, write
\[
\widetilde S_a=\kappa\log\sum_\tau e^{Q_{\tau a}/\kappa},\qquad
\widetilde V=\lambda\kappa\log\sum_a e^{\widetilde S_a/(\lambda\kappa)}.
\]
This is a valid sequential construction, but it is **not** a nested-GEV with fertility nests when \(\lambda<1\): that would require an invalid within/outer scale ratio \(1/\lambda>1\). Call it a “fertility-first sequential-logit control,” not a reversed nested logit. The document correctly warns it changes the random-utility law ([lines 85–89](file:///Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/docs/model/e5f_joint_nested_experiment.md#L85)).

3. **GEV law and primary probabilities are otherwise correct.** The CDF in lines 43–50 has marginal Gumbel scale \(\kappa\) and valid restriction \(0<\lambda\le1\); differentiating
\[
V=\kappa\log\sum_\tau\left(\sum_a e^{Q_{\tau a}/(\lambda\kappa)}\right)^\lambda
\]
with respect to \(Q_{\tau a}\) gives exactly the product probability in lines 57–60. At \(\lambda=1\), it reduces to flat four-way logit. Clarify that \(V\) is the *mean-zero-shock* expected maximum; under the stated uncentered CDF it is \(\gamma\kappa+V\) ([lines 52–60](file:///Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/docs/model/e5f_joint_nested_experiment.md#L52)).

4. **Feasibility needs a Boolean branch rule, not only dead values.** The intended rule is coherent ([lines 39–41](file:///Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/docs/model/e5f_joint_nested_experiment.md#L39)), but code must exclude an attempt if *any positive-probability* outcome is infeasible, and retain it for zero-probability infeasible branches. Do endpoint branching at \(\pi=0,1\); never evaluate \(0\times(-\infty)\) or let the current solver’s finite dead sentinel enter \(Q_{\tau1}\). Current tenure kernels use \(-10^{10}\) sentinels ([kernels lines 615, 651–656](file:///Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/intergen_eqscale_seq_optimized/kernels.py#L615)).

Future uncertainty does not break the four-plan interpretation if \(q(\tau,x)\) is the conditional expected continuation value and the new GEV shocks are current-date only. It does make this strictly a one-date perturbation, since frozen continuation retains the baseline future shock process, as the experiment appropriately states ([lines 9–14](file:///Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/docs/model/e5f_joint_nested_experiment.md#L9)).

No new joint-choice implementation or tests were present under `code/model`; I did not run model solves or inspect protected manuscript files beyond the requested read-only excerpt.