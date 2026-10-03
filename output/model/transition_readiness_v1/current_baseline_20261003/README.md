# Current-baseline fertility-shock estimation — paused before implementation

Tommaso paused this work on October 3, 2026, to resume at home. No transition
job was launched, no transition source was changed and no transition solve ran.
The separate ten-chain calibration array 19127370 continues with six-hour limits.

The author clarified the requested exercise: estimate one permanent 2007 change
in the fertility-preference parameter `psi_child` to match the retained
2020–2023 average fertility target, 1.64575. Earlier three four-year windows are
validation rows. A no-shock check is preliminary validation, not the deliverable.

The intended initial state is the verified current-parameter, one-birth Estate-A
case at `output/model/experiments/birth_count_choice/estate_a_v1/single/cases/20261003T212605039706Z_a739edc3`.
It retains post-interest transaction timing, soft financing and the current
chain-13 parameter vector, with estates net of the 6% selling cost in utility
and death-flow accounting, and no extra interest on post-saving financial wealth.
Pending recalibration candidates are not adopted initial conditions.

The maintained `transition_readiness/floor_runtime.py` constructor still loads
the older utility-floor/small-credit context. The legacy
`run_e5f_current_transition.py` additionally requires experimental natural
solvency. Neither can be silently used for this baseline. The existing
`one_shock_floor.py` controller supplies the scalar fitting loop, pinned targets,
checkpoints, budgets, two-horizon checks and strict numerical gates.

Next: verify dated backward and forward bindings against the current model,
then implement an isolated current-context runtime/controller adapter in the
birth-count experiment. Dynamic post-interest timing and net-estate accounting
are not yet verified. Preserve both birth-entry queues, inherited population,
fixed baseline housing-supply coefficient, fixed payroll tax/endogenous pension,
and the retained provisional estate-accounting limitations. Do not substitute
the historical fixed-stock path without reconciling the closure.

A preliminary design considered a six-hour, one-core cluster job with 96 GiB,
a 2-GiB exact-policy cache and two horizons of 24/32 four-year dates; this is
unreviewed and not launch-ready. Estimate its solve count and runtime, validate
the exact loop and current-baseline no-shock mapping, preserve all acceptance
checks, and label finite-horizon results diagnostic. No old Jacobian/checkpoint
can be reused as current evidence without verification.

Read-only grounding is in `work/inspect_final.md`. The interrupted implementation
agent `current_transition_engine` completed context loading only and changed no
files. Resume from this specification, avoiding a new broad historical audit.
