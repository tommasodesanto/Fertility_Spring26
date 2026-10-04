# Bounded read-only memory review

Verified 2026-10-04 by transition_contract_review. No model solve, source edit, or allocation profile was run. Sizes below were read from NPY headers in the selected saved native archive, not estimated from process RSS.

Reference: `output/model/experiments/birth_count_choice/estate_a_v1/single/cases/20261003T212605039706Z_a739edc3`.

## Exact saved array payload

- Core state array: `(120, 6, 1, 17, 9, 4, 4)`, float64, 14,100,480 bytes = 13.447265625 MiB.
- Each birth-count action/realized probability array has the same dimensions plus a four-outcome axis: 53.7890625 MiB each. The four-outcome allocation exists even with the birth cap set to one; the cap does not compress that axis.
- Tenure probability array: core state dimensions plus six outcomes, float32, 40.341796875 MiB.
- All saved parameter-object array payloads: 201.7224884033203 MiB; all saved solution array payloads: 358.6935272216797 MiB. These sums do not count Python object overhead or transient allocations.

## Retained copies and lifetime

`code/model/experiments/transition_readiness/floor_runtime.py:599-602` captures a copied pre-decision distribution and a deep-copied full dated parameter object for every date. Thus a current dated-state capture costs approximately 215.170 MiB/date if its arrays retain the saved shapes. This is inherited behavior, not a newly introduced horizon list. The count-specific parameter arrays add 107.578 MiB for the two count probability arrays and 40.342 MiB for three count distributions relative to an otherwise identical parameter object with those fields removed. That is an equal-state-space incremental comparison, not a measurement of the historical binary model's total memory or dimensions.

`code/model/tools/run_e5f_perfect_foresight_transition.py:594` retains H+1 value arrays, adding `(H+1)*13.447 MiB`. The forward/backward loops also have transient parameter copies. The static payload estimate for captured states plus values is approximately 5.37 GiB at H=24 and 7.16 GiB at H=32, or 12.53 GiB when both paths overlap. It excludes the baseline, policy cache, endpoint packets and transient solver allocations. This is a payload estimate, not a peak-RSS prediction.

`code/model/experiments/transition_readiness/one_shock_floor.py:809-818` uses one `latest` native mapping, overwritten after the next mapping completes. It does not keep every root evaluation. During construction of the replacement, the previous native mapping remains alive. The root evaluator returns only gates and residual vectors. `pinned_tools/e5f_social_security_root.py:124-138` passes those to the generic root; `pinned_tools/e5f_ssj_scaled_step_root.py:98-123` retains scalar/vector ledger records and current/best points, not full native policies.

`one_shock_floor.py:979-1016` retains the previous horizon until the next horizon can be compared. `:1038` retains only the latest completed candidate reply in `self.last`; it can overlap with the next candidate's work. `:744-784` intentionally memoizes one stationary endpoint packet per distinct successful psi, so endpoint storage grows with accepted distinct candidates, bounded by the run's candidate budget. `:667-676` warm snapshots contain paths and Jacobians, not full policy arrays.

The exact policy cache is scoped to one mapping and bounded by the current adapter's 2 GiB cap, compared with the inherited class's 64 GiB cap. No evidence of an unbounded full-path/root-evaluation leak was found. Capturing full parameter objects at every date is an avoidable retained-copy cost for the current date-4 comparison/export interface, but changing that implementation requires separate scoped verification; no change to tonight's running pinned source is proposed.

## Interpretation and limits

The smaller current state space does not by itself imply lower peak memory because duplicated arrays and overlapping path lifetimes contribute materially. The observed six-date smoke RSS of 13.93 GiB cannot be decomposed exactly using this static inspection. No OOM was observed or established: the reported v6 failure was the stationary endpoint's 16-evaluation convergence budget. The current seed-map observation reported by the lead (about 4 GiB RSS and about 100% CPU per process) is compatible with stage-dependent memory, but is not an allocation profile. These findings do not block the running launch.
