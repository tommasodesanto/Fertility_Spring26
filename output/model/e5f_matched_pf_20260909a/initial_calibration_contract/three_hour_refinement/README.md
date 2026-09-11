# Three-hour initial-calibration refinement

Author authorized useful cluster work during a three-hour teaching block on
September11. Job17375937 uses24 CPUs and192GiB, with a three-hour hard Slurm cap.
It runs numerical work without an active AI monitoring loop. No economic model,
empirical target, weight, parameter bound or numerical gate is changed.

This is the pre-2007 stationary working calibration, with new utility, sequential
fertility choices, payroll tax0.179, actual-budget pensions, housing elasticity
0.63 and fertility normalized separately to2.1. It is not a static2023 calibration,
a fitted historical transition or a policy benchmark.

## Search design and budget

The complete objective retains12 scored rows plus the separate2.1 normalization,
with nine structural coordinates. Fixed canonical objective:
`c0e266d3a0d430343c469d780d1aedb45fa87f8763c9c938889e0c37daa31de2`.
All641 source files and the original seed checkpoint are checked by the unchanged
wrapper before every model call. The economic snapshot remains70abd4a8. Its
exact solve/observe/score smoke17370427 passed two identical repetitions in979.51
seconds. The new orchestration has6 passing pure tests, including parallel
completion, failure preservation, transformed bounds and objective residuals.

1. Evaluate24 previously prepared joint proposals across all nine coordinates.
2. If all cases pass, take18 central derivative probes around the best completed
   candidate and propose12 jointly adjusted points from that fresh Jacobian.
3. Repeat step2 once if sufficient time remains.
4. Freeze selection and run two exact repetitions, including all original17plots.

Maximum84 search evaluations plus2 selected repetitions: up to688 stationary
GE calls, about344 at the observed four-solve normalization cost. Each case keeps
its existing1800-second solve and2100-second wrapper caps. With24 independent
workers, approximate duration50–100minutes excluding queue, with a3h hard cap.
Search has7800seconds; the last3000seconds are reserved for verification/reporting.
No stage starts without2100seconds left in its search window. Failed cases are
preserved and stop subsequent adaptive stages for review; there are no silent
retries, automatic target demotions or production promotion. Finishing before
three hours is allowed. A failed batch can still leave a useful better candidate.

## Saved evidence

Remote base:
`/scratch/td2248/projects/Fertility_Spring26_recent_parent_probe_70abd4a8/batches/three_hour_refinement_20260911`.

`results/heartbeat.json` updates every30seconds; per-solve progress remains in
`results/cases/<case>/evaluation/raw/heartbeat.json`. Each completed case updates
`results/latest_completed.json`, `best_so_far.json` and `cases.json`.
The final `summary.json` reports whether selected repetitions passed. All13
fit rows and17parameter rows are in `selected_target_fit.csv` and
`selected_parameters.csv`; all17 original plots are in
`selected_standard_diagnostics/`. Large model checkpoints remain on Torch.
Each case retains contracts, raw gates, source receipts, full fit tables and
parameter tables. Finite-difference matrices and proposed joint updates are saved.

The leading inherited candidate is retained even if new evaluations all worsen
fit. This working criterion still contains the documented resident-child observer
approximation and synthetic weighting scales; minimizing it does not certify
identification, fit quality, horizons or policies. No unattended AI heartbeat
has been reactivated, avoiding continuing token consumption while computation runs.
