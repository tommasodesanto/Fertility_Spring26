# September 28 morning fit comparisons

Supplemental fit figures for the authenticated selected overnight candidate
`initial_0506_block`. These do not replace the stable 17 standard diagnostics.
`build.py` loads the frozen saved checkpoint on Torch, verifies its receipt hash,
and imports the pinned runtime solely to deserialize and read income primitives.
It performs no equilibrium solve or calibration. One CPU, 32 GiB, 540-second
process cap and ten-minute Slurm allocation; job 18713061.

Outputs: `lifecycle_fit.png`, `targeted_fits.png`, `untargeted_fits.png`, combined
`fit_plots.pdf`, `model_profiles.csv`, and `qa.json`.

Empirical age profiles reuse only non-Model rows from
`../../daytime_calibration_20260927/lifecycle_dashboard/lifecycle_comparison.csv`;
its `provenance.json` documents samples and source verification. The historical
model rows are excluded. CPS June 2004/2006 profiles are capped children ever
born; complete four-year age bins only. Exact age 25 is a separate calibration
moment with within-period interpolation, not the 22–25 bin midpoint. Model
children profiles average pre/post birth counts and cap the last state at three;
completed fertility normalization uses a different top-bin weight.

ACS 2005/2006 all-structure household-head profiles use capped-nine rooms and
actual ownership, which differ from the aggregate AHS room and DUE ownership
calibration samples. PSID 2005/2007 net worth uses each sample's mean working-age
annual earnings denominator. These cross-sectional age profiles are descriptive
untargeted checks, not cohort trajectories or newly adopted targets. No
confidence intervals have been constructed. Connecting lines are for display.

The plotting job asserts exact replay of age-25 fertility and aggregate
wealth/earnings against the current selected target table. Original model,
parameters, weights, targets, samples, and scientific gates are unchanged.
