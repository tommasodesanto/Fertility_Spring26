# Income time aggregation diagnostic

Mode: `full`; batches: `20`; households per batch: `20000`; annual years: `120`.

The annual process is simulated with stationary persistent log variance and iid transitory log variance from the retained no-fixed candidate. Four-year levels are arithmetic averages of four annual mean-one levels. The 15-state and continuous endpoint objects are reported as approximations for comparison; neither is labeled an exact block-average process.

For annual mean-one levels, the exact covariance is `exp(Vp * rho^|d| + Ve * 1[d=0]) - 1`; block covariance is the average of its 16 annual pair covariances. Block log moments are reported from simulation because the log of an arithmetic average has no matching closed form here.

Receipt status: **completed**. Mean-one checks: **PASS**; exact level-moment checks: **PASS** under 6 batch standard errors + 0.001 absolute floor.

See `receipt.json` for complete moments, Monte Carlo SEs, source hashes, and validity flags. No raw panels are retained.

Supplemental plots: `level_covariance_comparison.png, block_average_quantiles.png`.
