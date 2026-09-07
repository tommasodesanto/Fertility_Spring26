Static audit complete; no source edits or model runs.

**Proven findings**

- The renter objective uses effective resources
  \[
  R^* = R_g b+y+\min\{\bar g,\max(\bar g-[R_g\max(b,0)+y],0)\},
  \]
  where \(\bar g=SD.g_b(n,cs)\). The transfer is debt-blind and must be included in reporting. See [solver.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/intergen_eqscale_seq_optimized/solver.py:2519) and [kernels.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/intergen_eqscale_seq_optimized/kernels.py:895).

- Let \(c_b=SD.c_b(n,cs)\), \(h_b=SD.h_b(n,cs)\), \(\alpha_v=SD.\alpha(n,cs)\), rent \(r_i\), and chosen saving \(b'\). The exact renter surplus is
  \[
  S=R^*-c_b-r_i h_b-b'.
  \]
  Feasibility in `eval_renter_scalar` is \(S>10^{-10}\); infeasible points receive \(-10^{10}\). [kernels.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/intergen_eqscale_seq_optimized/kernels.py:55)

- For a feasible uncapped solution,
  \[
  c=c_b+\alpha_v S,\qquad
  h=h_b+\frac{(1-\alpha_v)S}{r_i}.
  \]
  The cap applies iff
  \[
  h_b+\frac{(1-\alpha_v)S}{r_i}>h_R^{\max}
  \quad\Longleftrightarrow\quad
  S>\frac{r_i(h_R^{\max}-h_b)}{1-\alpha_v}.
  \]
  Then the objective uses
  \[
  c=c_b+\max\{R^*-c_b-r_i h_R^{\max}-b',\,10^{-10}\},
  \quad
  h=h_b+\max\{h_R^{\max}-h_b,\,10^{-10}\}.
  \]
  This is exactly the scalar/exhaustive objective, including its numerical \(10^{-10}\) guards. [kernels.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/intergen_eqscale_seq_optimized/kernels.py:59), [kernels.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/intergen_eqscale_seq_optimized/kernels.py:761), [kernels.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/intergen_eqscale_seq_optimized/kernels.py:966).

- Equivalence scales affect the reporting allocation through \(\alpha_v\), not by dividing \(c\) or \(h\) by the scale. `escale` multiplies utility only. Under the active eqscale specification, \(c_b=0\); \(h_b=0\) unless the child-room floor is active, while \(\alpha_v\) falls with at-home children unless that floor is active. [solver.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/intergen_eqscale_seq_optimized/solver.py:2260), [kernels.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/intergen_eqscale_seq_optimized/kernels.py:878).

- The observed defect is real: before the pending diff, feasible renter policies report \(c_b+\max(c_t,c_{\min})\) and \(h_b+\max(h_t,.01)\), although neither \(.04=c_{\min}\) nor \(.01\) is in the objective. [kernels.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/intergen_eqscale_seq_optimized/kernels.py:973)

**Review of the available lead diff**

- The diff is mathematically correct. For `exhaustive_saving and v_best > -1e9`, it replaces only those output floors with the precise \(c_t,h_t\) used by the objective, in both cap regimes. [kernels.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/intergen_eqscale_seq_optimized/kernels.py:978)

- Its eligibility mask is sufficient: `v_best > -1e9` implies a feasible scalar evaluation, hence \(S>10^{-10}\). The preceding infeasible branch remains unchanged. This prevents rewriting dead branches.

- Default-off sequential behavior remains byte-for-byte unchanged: the added assignments are gated by `exhaustive_saving`; the only active Markov caller supplies it as `int(joint_active)`, while the legacy core omits the optional argument and thus uses zero. [solver.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/intergen_eqscale_seq_optimized/solver.py:2528), [solver.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/intergen_eqscale_seq_optimized/solver.py:3028).

- Edge case: `h_b\ge h_R^{\max}` is ruled out for the active child-room-floor configuration by parameter validation, so the odd numerical cap guard is not live for the experimental target domain. [parameters.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/intergen_eqscale_seq_optimized/parameters.py:534)

**Fresh-history implication**

This cannot be treated as reporting-only for equilibrium outputs. At fixed prices it does not alter \(V\), \(b'\), or joint choice probabilities, but corrected renter \(h_R\) directly enters rental demand, housing residuals, policy statistics, event-study housing, and property-tax revenue. A fresh equilibrium/history evaluation is therefore required for any reported result. [solver.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/intergen_eqscale_seq_optimized/solver.py:7434)

No conjectural calibration conclusion.