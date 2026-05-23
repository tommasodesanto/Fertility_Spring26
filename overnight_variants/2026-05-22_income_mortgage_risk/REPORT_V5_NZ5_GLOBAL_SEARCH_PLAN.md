# V5 Nz=5 Global Search Plan

This plan launches a true global randomized search for the isolated V5
benchmark-normalized outside-option HANK-\(z\) prototype.

## Scope

- Model: V5 outside-option closure with benchmark normalization inside the GE
  loop.
- Grid: `Nb=30`, `Nz=5`, \(\rho_z=0.95\), unconditional \(\sigma_z=0.35\).
- Objective: full 19-moment weighted SMM loss, plus penalties for non-converged
  GE or invalid outside-option scale accounting.
- Workers: `8` parallel Python worker processes.
- Thread cap: `OMP_NUM_THREADS=1`, `MKL_NUM_THREADS=1`,
  `OPENBLAS_NUM_THREADS=1`, `NUMBA_NUM_THREADS=1`.
- Time budget: `28800` seconds, approximately 8 hours.

At roughly `100-300` seconds per candidate under parallel contention, the
expected run size is about `800-2300` completed evaluations. The wide interval
is deliberate because first-time compilation and eight-process contention can
move realized solve times substantially.

## Parameters

The global search samples 19 parameters:

`beta`, `b_entry_fixed`, `psi_child`, `h_bar_jump`, `h_bar_n`, `c_bar_n`,
`kappa_fert`, `chi`, `kappa_loc`, `mu_move`, `theta0`, `theta_n`, `h_bar_0`,
`E_C`, `r_bar_C`, `alpha_cons`, `phi`, `hR_max`, and `h_own_max`.

The seed bank includes the current V5 `Nz=5` baseline and the best directional
probes, but after the seed bank the proposal is mostly global: the default
`global_prob=0.75` draws uniformly over the full bounded hyperrectangle, while
the remaining draws perturb the current best.

## Checkpointing

Each run writes to `global_search_v5_nz5/<run_tag>/`.

Files written continuously:

- `evaluations.jsonl`: full record for every completed candidate.
- `evaluations.csv`: compact candidate/moment table.
- `latest.json`: latest completed candidate.
- `best.json`: current best objective candidate.
- `best_summary.md`: readable current best parameter and moment table.
- `status.json`: submitted/completed/active counts and best loss.
- `heartbeat.txt`: simple text heartbeat for shell monitoring.
- `driver.log`: stdout/stderr from the background driver process.
- `driver.pid`: background process id.

## Smoke Test

The exact parallel/checkpointing loop was smoke-tested with two workers and two
small-grid evaluations. It produced `evaluations.jsonl`, `evaluations.csv`,
`best.json`, `best_summary.md`, `latest.json`, `status.json`, and
`heartbeat.txt`, then exited normally. The smoke penalties were expected because
the smoke capped the GE loop at two iterations.

## Launch Command

```bash
cd /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/overnight_variants/2026-05-22_income_mortgage_risk
RUN_TAG=v5_nz5_global_$(date +%Y%m%d_%H%M%S)
mkdir -p "global_search_v5_nz5/${RUN_TAG}"
OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMBA_NUM_THREADS=1 \
nohup caffeinate -dimsu /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/.venv/bin/python \
  run_v5_nz5_global_search.py \
  --run-tag "${RUN_TAG}" \
  --results-dir global_search_v5_nz5 \
  --workers 8 \
  --budget-sec 28800 \
  --max-evals 1000000 \
  --seed 20260523 \
  --nb 30 \
  --nz 5 \
  --rho-z 0.95 \
  --sigma-z 0.35 \
  --kappa-entry 1000000 \
  --max-iter-eq 35 \
  --tol-eq 5e-4 \
  --global-prob 0.75 \
  --log-every 10 \
  > "global_search_v5_nz5/${RUN_TAG}/driver.log" 2>&1 &
echo $! > "global_search_v5_nz5/${RUN_TAG}/driver.pid"
```

If `caffeinate` is unavailable, remove `caffeinate -dimsu` from the command.
