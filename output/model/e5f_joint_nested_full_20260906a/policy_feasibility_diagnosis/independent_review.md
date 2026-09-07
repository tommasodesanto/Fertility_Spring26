Static diagnosis: the proposed unsupported-endpoint mechanism is real in the joint path. It is not yet proof that it generated the two reported cells; the lead’s saved-state replay should establish that final causal link.

- [`kernels.py`](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/intergen_eqscale_seq_optimized/kernels.py:573) linearly interpolates a finite dead sentinel, \(-10^{10}\), in `_interp_with_clip`. `tenure_choice_kernel` uses it for an owner selling to renter at [631](</Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/intergen_eqscale_seq_optimized/kernels.py:631>), renter buying at [646](</Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/intergen_eqscale_seq_optimized/kernels.py:646>) / [654](</Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/intergen_eqscale_seq_optimized/kernels.py:654>), and owner switching products at [659](</Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/intergen_eqscale_seq_optimized/kernels.py:659>). No endpoint-support test occurs.

- Thus, with adjacent conditional values \(V_k\le -10^9\) and \(V_{k+1}>-10^9\), positive weights \(1-w,w\), the code admits
  \[
  \tilde V=(1-w)V_k+wV_{k+1}.
  \]
  Since the sentinel is finite, \(\tilde V>-10^9\) whenever the dead-endpoint weight is sufficiently small (roughly \(<0.1\), ignoring the finite live value). It can therefore defeat the `q > -1e9` alive test in [`joint_nested.py:73`](</Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/intergen_eqscale_seq_optimized/joint_nested.py:73>).

- The joint Bellman then records that product in `joint.products` ([74](</Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/intergen_eqscale_seq_optimized/joint_nested.py:74>), [2636–2643](</Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/intergen_eqscale_seq_optimized/solver.py:2636>)). `factor_age` accepts it because it checks only product-code validity, not conditional-value support ([143–153](</Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/intergen_eqscale_seq_optimized/joint_nested.py:143>)).

- The forward transaction map then scatters the selected mass across precisely those wealth endpoints ([3523–3525](</Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/intergen_eqscale_seq_optimized/solver.py:3523>); scatter at [3600–3605](</Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/intergen_eqscale_seq_optimized/solver.py:3600>)). This is the discrete-grid leak: a selected branch uses a conditional policy that is dead at one positive-mass destination node.

The stated gate limitation is also confirmed. The transition driver gates pre-fertility mass against aggregate `policy.V` ([366–396](</Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/run_dynamic_population_transition.py:366>)), while the problematic object is the destination-tenure conditional \(V_d(\cdot,\tau)\). A destination renter may be conditionally dead yet have aggregate joint value alive because another tenure plan is feasible. Hence the gate cannot enforce support for the actually selected renter branch. This is not fixed by `_censor_entry_dead_mass`, which uses the same aggregate value criterion ([3938–3964](</Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/intergen_eqscale_seq_optimized/solver.py:3938>)).

The reported \(c=0.04\) is consistent with this route: the compiled renter evaluator assigns the reporting floor on `surplus <= 1e-10` while retaining `-1e10` value ([540–547](</Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/intergen_eqscale_seq_optimized/kernels.py:540>). That is a dead conditional renter output, not evidence that the already-repaired feasible-renter reporting rule is wrong.

Appropriate narrow repair if replay confirms it:

\[
\text{admissible}_{\tau}(x)
\iff
\bigl[(1-w)>0\Rightarrow V_{k,\tau}>-10^9\bigr]
\land
\bigl[w>0\Rightarrow V_{k+1,\tau}>-10^9\bigr].
\]

If false, set that branch’s value to \(-10^{10}\) before the tenure argmax / joint nest. This is a numerical support restriction for interpolation, not a model primitive or a post-hoc mass deletion. Make it an explicit default-off `strict_interpolated_support` argument used only by `joint_nested.bellman_block`’s call to `tenure_choice_kernel`; the sequential call at [`solver.py:3260`](</Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/intergen_eqscale_seq_optimized/solver.py:3260>) remains byte-identical.

Trace method for the two cells:

1. At the predecessor cell selected into the renter branch, record origin \((b,to,j,z,n,c)\), product \(tn=0\), post-transaction \(x=b+(1-\psi)pH_{to}\), interpolation \((k,w)\), both `Vd[k,0,...]`, `Vd[k+1,0,...]`, aggregate joint `V`, and selected joint probability.
2. Verify \(0<w<1\), one renter endpoint is \(\le-10^9\), aggregate joint \(V>-10^9\), and the other endpoint is live.
3. Reconstruct the exact forwarded contributions \(m(1-w)\) and \(mw\) from the tenure map. One must equal the reported renter endpoint mass, up to accumulation roundoff.
4. Confirm the dead renter evaluator has `surplus <= 1e-10` and produces its default output floor. This distinguishes the leak from a feasible-policy reporting defect.

Pitfalls for implementation/testing:

- Treat zero weights as non-supporting: clipped lower/upper endpoints and exact-grid \(w=0\) or \(w=1\) must inspect only the endpoint receiving mass.
- Apply the condition after constructing the same branch wealth as the current code, including sale proceeds, grants, and the \( (1-\phi)pH \) down-payment / \(-\phi pH\) borrowing constraints. Do not change those constraints.
- Cover every interpolated tenure branch, including owner \(\rightarrow\) renter, not only renter \(\rightarrow\) owner; direct same-tenure values need no endpoint test.
- Do not rely on argmax ties: all-dead menus intentionally retain default product code but receive zero joint probability after the `q > -1e9` filter.
- The stochastic `tenure_logit_kernel` duplicates the same interpolation pattern ([703](</Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/intergen_eqscale_seq_optimized/kernels.py:703>) onward), but the joint Bellman currently calls deterministic `tenure_choice_kernel`. Keep any change scoped to the confirmed joint route unless that separate path is deliberately audited.

No files, jobs, or model runs were changed.