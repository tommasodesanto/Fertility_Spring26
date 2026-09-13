# Claude handoff: sequence-space transition root

You are continuing a narrowly scoped SSJ investigation in
`/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26`. The goal
is a usable accelerated transition root with honest benchmarks. Use the existing Claude session and Torch for numerical work; do not call
paid external model APIs or touch the independent four-shock transition arm. Do not edit the protected
manuscript subtree `latex/JMP_DS_draft/`.

Start with read-only inspection of `memory/AGENT_MEMORY.md`, the latest available
`memory/daily/YYYY-MM-DD.md`, `CALIBRATION_STATUS.md`, and then the
following SSJ artifacts:

- `code/model/tools/e5f_sequence_space_prototype.py`
- `code/model/tools/test_e5f_sequence_space_prototype.py`
- `code/cluster/check_e5f_sequence_space_native.py`
- `docs/model/e5f_sequence_space_prototype.md`
- `output/model/e5f_sequence_space_prototype_20260913/native_smoke/retry2/native_smoke.json`
- `output/model/e5f_sequence_space_prototype_20260913/native_smoke/retry2/lead_verification.json`

The native evidence is six evaluations over two dates in 117.57 seconds. The
baseline scaled residual is $7.6768292368\times10^{-5}$, below the unchanged
$2\times10^{-4}$ gate. State replay and residual reconstruction are exact;
the stationary drift checks pass. Halving the log-price step changes the
directional derivative by 0.1578425% in relative sup norm, with a maximum
substantial-component difference of 0.300674%. This is a derivative smoke,
not toolkit integration or a demonstrated speedup.

The mathematical unknowns are `(log_house_price, pension, rebate)` at each
date. The native residual remains housing relative imbalance plus the existing
PAYGO and equal property-tax-rebate relative imbalances, with fiscal entries
scaled by 200 exactly as documented. Preserve the fixed native household
kernels, population law, original birth queue (births/2.1), no immigration,
equal property-tax rebate, PAYGO, all feasibility/mass/budget gates, and the
native level-valued carried state `(g_pre, scheduled_entries,
scheduled_raw_entries)`. Never substitute the terminal stationary state.

Use the frozen checker bindings in
`code/cluster/check_e5f_sequence_space_native.py`: `c.rebated`,
`c.queue.queue_path`, `rebated.InheritedState(2007, c.old.initial_state)`,
`NS(..., asset_price=q)`, and `rebated.dated_residual` followed by
`rebated.stack_dated_residuals`. The authoritative remote evidence directory
is `torch:/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/e5f_sequence_space_native_20260913_retry2/output`.

Proceed in this order: inspect the two-date evidence; then build and test a
bounded Jacobian/root experiment locally in coordinates, using Torch compute
where it actually reduces native mapping work. Do not launch a blind dense
100-period finite-difference sweep. Compare an identical baseline, tolerances,
and target closure against the existing root; report actual wall times,
residuals, state/gate checks, and derivative directions. If a change is
warranted, implement it in an isolated, clearly named prototype with focused
tests and a full snapshot plus profiling/checkpoint artifacts. Preserve all
production files. Give an initial assessment within 30 minutes, then smoke-test
and cap the first numerical prototype at two hours, with checkpoints and a clear
stop condition. No paid agent supervision loops or unbounded cluster jobs. The lead agent is handling the
four-shock launch now, so do not cancel or modify it.
