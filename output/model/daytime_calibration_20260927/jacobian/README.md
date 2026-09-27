# Local calibration Jacobian — September 27 daytime

This is a diagnostic of the overnight main winner `de_0093` (loss 42.281937),
not a new specification or a claim of optimizer convergence. It reruns the
same internally normalized objective at eighteen central-difference points.
The child-benefit level is normalized separately at every point, exactly as in
calibration. All source, targets, weights, numerical gates, grids and bounds
are retained.

- Live output: `run_v2/`; `heartbeat.json` every five seconds, and
  `checkpoint.json`, `latest_completed.json`, `best_so_far.json` every case.
- Frozen plan and contract: `tmp/daytime_calibration_20260927/jacobian/` from
  repository root. Plan records complete coordinates, physical point values,
  runner SHA, source-contract ancestry and current contract SHA.
- Controller: `code/model/tools/run_e5f_local_jacobian.py`.
- Budget: two local evaluators, eighteen cases, ninety minutes total,
  fifteen minutes per case. The first minus/plus H0 pair must both pass before
  any other coordinates are dispatched. Failure halts new dispatch; already
  running cases finish within their existing caps. No automatic retries or
  additional step-halving cases.
- Each case stores its full fourteen-row target fit and 31-row parameter table,
  normalized benefit, scientific receipt and saved stationary solution. The
  existing immutable evaluator validates these receipts and hashes.

The six positive coordinates H0, chi, first-birth fixed cost, both fertility
scales and theta0 use z = log(parameter), step 0.02. The three remaining
coordinates are beta_annual / 0.05, delta_alpha_jump / 0.25, and
child_benefit_curvature / 0.8, each with step 0.01. These definitions are local
sensitivity units, not new search bounds. Every plus/minus point lies within
its original parameter bounds.

On successful completion, `jacobian.csv` gives all 14 x 9 unweighted
central derivatives and second differences; `jacobian.json` includes raw and
sqrt(weight)-scaled matrices, singular values/vectors, the local loss gradient
and normalized child-benefit responses. The normalization row is reported but
has zero scoring weight. Singular values depend on these declared coordinate
units. Grid kinks, normalization tolerance and untested step size can limit
interpretation; do not present the numerical rank as global identification.

`run_v1/` is an empty background-start attempt: no evaluator or solution was
launched there. It is preserved. `run_v2/` is the actual foreground-managed run.
