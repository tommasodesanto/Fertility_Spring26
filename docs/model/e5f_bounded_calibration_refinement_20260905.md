# Bounded same-contract calibration refinement, 5 September 2026

At 15:27 UTC the author authorized continuing the proposed candidate verification
and a small joint refinement while he works on theory. This is a new authorization
after the completed overnight audit. It does not change the model, target system,
weights, parameter bounds, or external restrictions. The existing production
benchmark remains September 4 task_010 until the new results are reviewed.

## Design and economics

The model overpredicts mean occupied rooms and underpredicts rooms gained after a
first birth. Search all eleven existing parameter coordinates against the complete
twelve-row objective; report both housing moments and every other target for every
candidate. No target dropping, reweighting, activated young-ownership row, imposed
bequest restriction, or rental-menu change is part of this run.

The original 23-case diagnostic panel contains a better twelve-row candidate:
task_003, loss 29.385883689710358, compared with production 30.482966707698903.
This is the original radius-0.02 panel, **not** task_003 in the later radius-0.005
small-step panel. Re-run the original generator twice with the production
twelve-row source. Require exact reproduction of all twelve fit rows, estimated
parameters, preference normalization, and saved numeric historical path entries.
The two runs smoke-test the complete calibration, observer, collection, checkpoint,
and seventeen-graph path before any expansion.

Round 1: a complete 23-case coordinate panel around the verified better point,
with radius 0.005 in the existing bounded transformed coordinates. Eleven free
parameters and twelve informative rows are retained. Weak local identification
remains a limitation; neither a smaller objective nor a numerical derivative
establishes precise parameter identification.

Round 2: at most twelve jointly varied candidates, chosen by a documented direct
pattern-search rule from round 1. No inverse-Jacobian or ridge proposal is used;
the original ridge planner and its rejection remain unchanged. Predeclare the
joint directions and step radii before launching this round. Retain the best
observed valid point even if every new joint candidate worsens the objective.
Use at most two exact repeats of the final best point. Stop after round 2.

The round-2 rule, set before round 1 completes, starts at its lowest-loss point.
One family combines all individually improving coordinate moves; another combines
the three with the largest full-objective improvement. Each uses one-half, one,
and two times the original 0.005 coordinate move. Two further families increase
the first-child room intercept by 0.025, 0.05, or 0.10 rooms and reduce the
per-child floor by one third of that amount, either alone or together with the
successful non-floor coordinate moves. Since the coded floor is
\(\bar h(m)=\bar h_J+m\bar h_C\) for \(m>0\), this preserves the three-child
floor and raises the first-child floor by two thirds of the intercept change.
These are explicitly joint economic parameter probes, not estimates of a local
derivative. All existing parameter bounds are enforced; duplicate patterns are
omitted. No numerical Jacobian inverse is needed to evaluate these hypotheses.

## Budgets and gates

- Overall ceiling: 19:27 UTC (four hours from authorization), including reporting.
  New candidate launches stop at 18:27 UTC; reserve the final hour for collection,
  verification, report rendering, and source backup. Queue delays are counted.
- At most 39 full historical evaluations: two initial repeats, 23 coordinate
  cases, twelve joint cases, and two final repeats. One CPU and 8 GB per case,
  at most six concurrent cases, 30 minutes per case. Observed predecessor costs
  were approximately 7–11 minutes per case: 4.6–7.2 serial CPU-hours for all 39,
  approximately 1–2 hours of computation at six-way concurrency plus queueing.
  Each evaluation normalizes the old stationary state and clears five historical
  markets; the driver records actual backward-solve counts rather than treating
  an evaluation as one Bellman solve.
- Smoke passage is required before expanding. If a scientific fingerprint,
  empirical input, history, accounting, feasibility, probability, or population
  check fails, retain the failure and investigate before additional dependent work.
- Preserve the production market gate 2e-4, accounting/population gates 2e-10,
  stationary stock/flow measurement gate 2e-8, and all driver gates. The separate
  ridge planner's stricter stock/flow gate is neither called nor changed.
- A one-minute wall-clock heartbeat accompanies each solve. Every completed case
  writes its complete fit, parameters, terminal checkpoint, and seventeen standard
  figures. Keep a latest-completed and best-so-far summary across the run.
  Investigate a running case without a heartbeat or checkpoint for 30 minutes.
- Stop on negligible improvement after the declared rounds; do not keep searching
  for a preferred policy effect. No policy or long-run transition is included.

## Execution and provenance

Use the production-matching snapshot at
`/scratch/td2248/projects/Fertility_Spring26_independent_audit_20260905`,
scientific bundle `630ba20bca6a1b54eb4c46aca904c4a087afb8c808b9c7f4660d5fcd316a970e`,
target fingerprint `3726c17e62c8233ce62d5f4c95f44fd2cc2ea6cfa3d2492795461b4569300497`,
and account `torch_pr_570_general`. The adapter only copies the terminal observer
arguments after the original callback. It does not change any solver, measurement,
transition, or target file. Its own hash, existing diagnostic helper hash, input
hashes, and complete launch arguments are pinned separately in each stage plan.

Code: `code/model/tools/run_e5f_bounded_calibration_refinement.py` and
`code/cluster/submit_e5f_bounded_calibration_refinement.sh`.
Results: `output/model/e5f_morning_refinement_20260905a/`.
The complete report will distinguish the retained production benchmark, verified
new candidates, and outstanding age-measurement/identification decisions.
All author-owned theory files and other agents' working changes remain outside
this task.
