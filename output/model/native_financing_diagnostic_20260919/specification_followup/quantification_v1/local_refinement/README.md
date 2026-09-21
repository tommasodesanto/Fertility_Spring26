# Twenty-minute local earnings refinement

Author authorization, September 21: refine locally for twenty minutes with ten
worker processes while discussion continues. This keeps the existing diagnostic
15-state persistent-plus-transitory candidate, all floors, nine free structural
coordinates, twelve scored target rows plus separate fertility normalization,
and original bounds (annual beta at most 0.99). No policy or specification change.

Starting point: selected sensitivity-panel objective 326.9831988727637;
`anchor_score.json` preserves its complete enriched score and parameter receipt.
Initial normalization seed: 0.23949950404168222.

Plan: validate relocated immutable inputs, then start a 1,200-second global
budget. First perform a two-repetition local anchor smoke using the actual
scored evaluator. Stop if the smoke fails. After successful smoke, run up to
20 deterministic nearby joint proposals with at most ten simultaneous workers,
each with one numerical thread. At most 24 full objective evaluations including
two smoke evaluations and two selected verification evaluations; at most 192
nested stationary solves. Selection verification runs only within the same time
budget; if time is insufficient, label the best new point provisional.

Recent cluster points took several minutes each; local concurrency and initial
compilation may make them slower. Expect a small number of waves, not a completed
optimization. The local smoke will supply the first machine-specific timing.
The frozen September 14 reference remains untouched. Exact source inventory,
income and target fingerprints are required; path relocation is documented.
Latest and best receipts and heartbeats remain accessible during execution.
Incomplete or rejected evaluations are reported, not assigned fabricated losses.

Status: launched; see `launch.json`, `controller.log`, and `run/heartbeat.json`.
The two-repetition local smoke gates all proposals. Three mocked controller
checks passed, including a failing-smoke stop, and native zero-solve wrapper
preflight passed. All 641 scientific source files match their pins.

The local isolated Python runtime reuses the project environment packages.
`checkpoint_compat.py` maps NumPy 2 and Python 3.13 pickle module names to the
local equivalent classes so saved checkpoints can load; array data and solver
code are unchanged. The smoke requires exact local repetitions and exact
numeric fit versus the saved anchor; no scientific gate is relaxed.

Proposal steps are up to 5% and 2% of anchor magnitudes, except beta uses
absolute steps up to .001 and .0004, clipped to existing bounds. This is a
finite nearby search, not a claim of convergence. Runtime staging is recorded
in `local_refinement_plan.json` and `relocation_receipt.json`.

