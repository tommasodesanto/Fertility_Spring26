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

## Completed and independently verified

All eighteen cases completed successfully in 50.57 minutes, within the 90-minute
budget. The first H0 pair passed the full evaluator checks before the remaining
sixteen were dispatched. The independent collection reran source/target pin
checks and the original evaluator's success validation for all eighteen cases:
checkpoint and receipt hashes, scientific identity, normalization, fourteen
complete target rows, 31 parameter rows and every searched bound. There are
252 target rows and 558 parameter rows in the combined tables. The central
matrix, singular values and weighted loss gradient were independently recomputed.
See `run_v2/independent_verification.json`.

The best diagnostic perturbation lowers loss from 42.282 to 41.993 by reducing
the first-birth fixed cost from 0.660 to 0.647. Reducing the initial fertility
taste scale from 0.195 to 0.191 gives 41.994; raising curvature from 0.071 to
0.079 gives 42.099; raising the continuation fertility scale from 0.347 to
0.354 gives 42.225. These demonstrate remaining improvement at the anchor.
They are diagnostic points, not independently repeated new calibration winners.

With the explicit local coordinate scaling above, the weighted 14-by-9 matrix
has numerical rank nine; singular values range from 317.614 to 1.210, with
condition number 262.387. This is only finite-difference local sensitivity.
The least-sensitive joint direction is mostly child-benefit curvature and the
first-birth fixed cost (right-vector coefficients -0.856 and -0.473), with
smaller initial fertility-scale and bequest coefficients. The bequest parameter
has the smallest individual weighted column norm (1.750), versus 315.057 for
chi. Signs of singular vectors are arbitrary. None of this establishes global
identification or precise statistical standard errors.

**Important step-size limitation:** beta and chi are not reliable inputs to an
unrestricted linear Newton step at these finite steps. The weighted symmetric
second difference, divided by the weighted plus/minus change, is 1.046 for beta,
0.318 for H0 and 0.151 for chi. For chi, the direct loss difference gives a slope
of -553.877 while `2 J' residual` gives -20.850: nonlinear effects matter even
though the moment-level second-difference ratio looks more moderate. Both chi
perturbations worsen loss sharply (92.402 and 70.247). This is evidence to
inspect narrower steps/grid kinks, not permission to extrapolate the central
slope. No half-step solves were run. Other ratios are below 0.019, except the
housing-loading ratio 0.005, and their loss slopes generally agree well.

The local fertility directions move early fertility only modestly. Increasing
the continuation taste scale by 2% raises children born by age 25 from 0.528 to
0.530; reducing the initial taste scale by 2% raises it to 0.529. The target is
0.810. This is local evidence, not proof that the target is unattainable.
Normalization remains near 2.100, with residuals roughly 0.00007--0.00013. Tiny
effects, such as the approximately 0.000006 plus/minus change in early fertility
from theta0, should not receive economic interpretation before tighter-noise
and step-size checks. These tolerances are warnings, not formal error bounds
on every moment.

Full supporting files in `run_v2/`:
- `all_target_fits.csv`: every target, model value, gap, weight and contribution.
- `all_parameters.csv`: all parameters, estimates, bounds and restrictions.
- `jacobian.csv` / `jacobian.json`: central derivatives, scaled SVD and gradient.
- `step_diagnostics.json`: finite-step nonlinearity and direct loss slopes.
- `independent_verification.json`: authenticated checks for every case.

No further solve or model change was launched during collection.
