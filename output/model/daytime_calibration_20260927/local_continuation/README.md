# Small local calibration continuation

COMPLETE13:07:38 EDT. All8 evaluations passed, no failures/timeouts/unrun
points. Two final repeats and export reproduce all14/31 table rows byte exactly;
all17plots visually reviewed. Best41.111684 versus anchor41.992733 and
original42.281937. No active workers. See run_v1/lead_completion_review.json
and run_v1/selected_export for full evidence. This is a small improvement,
not an optimum or complete policy-shape certification.

Explicit September27 author authorization: run a little locally to see how calibration evolves, after the cluster job never started. This new bounded run does not revive or extend the cancelled cluster schedule.

Launched12:45:55 EDT, controller PID2811, two actual single-thread workers. Search cutoff13:03:47 EDT; absolute end13:14:47 EDT including final repetitions and export. Eight-minute per-objective cap,29-minute global cap. Memory at launch:48GiB physical,82% free; no other local model processes active.

Immutable `runner_v1.py` and `plan_v1.json` select four nearby proposals using the unchanged authenticated frozen evaluator. Two anchor smokes first repeat the promising Jacobian first-birth-cost decrease (loss41.993 versus original42.282), comparing every numeric field of all14targets and31parameters exactly. Four proposals further lower first-birth cost, lower initial fertility choice scale, increase child-benefit curvature, or combine the last two. They avoid the beta/H0/chi directions with large finite-step nonlinearity. Bounds verified before launch; no weights, targets, economic rules or grids changed. The child-benefit intercept is normalized by the unchanged calibration procedure, as in the original objective.

At most8 objectives:2smokes,4search points,2fresh best-point repeats with17standard plots. Based on6stationary solves and346seconds per objective, expect roughly24minutes at2workers; hardstop remains29minutes. No objective is promised to complete beyond its budget. Fatal errors stop new dispatch, numerical inadmissibility and owned timeouts remain separately recorded, and unrun points are not assigned fabricated losses.

Read `run_v1/heartbeat.json`, `checkpoint.json`, `latest_completed.json`, `best_so_far.json`, `complete.json` and eventual `selected_export/` for progress and full fit/parameter tables. Source and all reference contract/checkpoint/table hashes are pinned. Parent signal/deadline handling closes owned child processes. Four focused controller tests and authenticated no-solve preflight passed. Actual exact-loop smokes are part of the run and gate proposal dispatch.

This is a local neighborhood exploration, not proof of optimum or identification. Native transition integration is separate and remains indexed in `../credit_transition/README.md`.
