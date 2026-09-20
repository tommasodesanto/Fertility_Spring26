# Income grid-resolution diagnostic

This deterministic algebra audit compares 5/7/9/15/25 persistent states crossed with 3/5 iid-transitory quadrature nodes at the retained annual parameters. It does not run the household model, simulate new panels, estimate a process, or adopt a process. The 5x3 row reproduces the existing constructor and the prior full MC payload fingerprint.

Grid `level_cov_lag*` and `log_cov_lag*` are endpoint Markov moments. The MC reference is the exact annual four-year block-average experiment already collected; its tail quantiles are therefore not directly comparable to endpoint-grid quantiles. The 25x5 variance happens to sit near the block-average variance through opposing discretization errors; refinement converges toward the continuous endpoint target, not the block-average target. The second plot displays that distinction explicitly.

Receipt status: **completed**. See `receipt.json` and `grid_moments.csv`; plots are `grid_resolution_moments_quantiles.png, grid_resolution_tail_comparison.png`.
