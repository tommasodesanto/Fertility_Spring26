# Six-hour M5 continuation

This workflow asks whether additional computation improves the live M5
15-moment objective. It remains inside the isolated optimized package and does
not redirect any production import.

Two baselines remain explicit. The canonical production M5 loss is
`9.044422069071352`. Re-solving the same theta with the optimized strict-root
evaluator in the promotion panel gives `9.14507565057346`, because the solvers
accept slightly different nearby clearing prices. The collector reports
improvement against both; beating only the latter is not described as beating
canonical M5.

## Contract

- Model and targets: unchanged M5 income-disciplined specification.
- Identification: 14 free parameters, 15 hard moments, `theta_n=0` externally.
- Search evaluator: `J=17`, `Nb=120`, `max_iter_eq=10`, `tol_eq=1e-4`.
- Promotion evaluator: two fresh repeats at `max_iter_eq=40`, `tol_eq=2.5e-5`.
- Search: eight one-CPU chains, alternating transformed Nelder--Mead and
  coordinate-pattern search, with predeclared step and start perturbations.
- Budget per chain: 340 search minutes, 10 strict-repeat minutes, 2,000
  evaluations, and a six-hour Slurm wall clock.
- Progress: every candidate appends to `cases.jsonl` and refreshes
  `latest_completed_case.json`; improvements refresh `best_so_far.json`.
- Stop: wall-time budget, evaluation cap, or pattern-step convergence.

The Torch promotion panel measured 14.2 seconds per optimized full-GE
candidate. The design therefore targets roughly 1,400--1,500 candidates per
chain, or approximately 11,000--12,000 candidates across eight chains, before
overhead and infeasible proposals.

## Exact-loop smoke

The smoke runs two array tasks and completes the full 14-dimensional initial
loop (15 evaluations) under a smaller numerical grid:

```bash
OPTIMIZED_M5_SMOKE=1 OPTIMIZED_M5_RUN_TAG=intergen_optimized_m5_continuation_smoke_20260719 \
  sbatch --array=1-2%2 --time=00:20:00 submit_six_hour_calibration.sh
```

Both smoke tasks must write `metadata.json`, `cases.jsonl`,
`latest_completed_case.json`, `best_so_far.json`, and `summary.json`, with no
unexpected program errors.

## Production-sized launch

```bash
OPTIMIZED_M5_RUN_TAG=intergen_optimized_m5_continuation_20260719 \
  sbatch submit_six_hour_calibration.sh
```

After obtaining the array job ID, submit the collector with
`--dependency=afterok:<array_job_id>`. The collector selects only chains with
two exact strict repeats and writes the complete target-fit, parameter/bounds,
chain, and result tables under the run's `report/` directory. Loose search
losses are never eligible for comparison to the canonical M5 loss.
