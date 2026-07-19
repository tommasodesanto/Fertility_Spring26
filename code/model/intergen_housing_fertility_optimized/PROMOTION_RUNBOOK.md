# Promotion battery runbook

This battery is isolated from production imports and writes only under
`output/model/intergen_optimized_promotion_20260719/`. It is a numerical
promotion exercise, not a new calibration or model specification.

## Declared workload

| Arm | Tasks | Work per task | Expected wall time |
|---|---:|---:|---:|
| Exact-bound parity | 29 | production + optimized fixed-price solve | under 5 min |
| Full demand map | 17 | production + optimized fixed-price solve | under 5 min |
| Actual-root cases | 6 | production + optimized strict GE solve | under 10 min |
| Strict repeats | 2 | one optimized strict GE solve | under 5 min |
| Diagnostics | 1 | two fixed-price solves plus full packets | under 10 min |
| Calibration throughput | 2 | ten sequential GE objective evaluations | under 30 min |

Every array task writes one atomic completion record. The calibration tasks
rewrite their checkpoint after every candidate. The wall-time cap is 45
minutes per task; no task should be silently restarted if that cap is reached.

## Local exact-loop preflight

From `code/model`, use `Nb=30` while preserving all loop branches:

```bash
PYTHONPATH=$PWD .venv/bin/python -m intergen_housing_fertility_optimized.promotion_worker \
  bound-parity --index 1 --output ../../output/model/intergen_optimized_promotion_20260719_smoke/bounds/task_1.json --smoke

PYTHONPATH=$PWD .venv/bin/python -m intergen_housing_fertility_optimized.promotion_worker \
  demand-point --index 10 --output ../../output/model/intergen_optimized_promotion_20260719_smoke/demand/task_10.json --smoke

PYTHONPATH=$PWD .venv/bin/python -m intergen_housing_fertility_optimized.promotion_worker \
  root-case --index 1 --output ../../output/model/intergen_optimized_promotion_20260719_smoke/roots/task_1.json
```

## Torch launch

First run the normal authentication/queue probe:

```bash
code/cluster/torch.sh status
```

From the remote project root, submit the exact loop smoke before the full
battery:

```bash
cd code/model/intergen_housing_fertility_optimized/cluster
PROMOTION_TASK=bound-parity PROMOTION_SMOKE=1 sbatch --array=1-2 submit_promotion_battery.sh
PROMOTION_TASK=demand-point PROMOTION_SMOKE=1 sbatch --array=9-11 submit_promotion_battery.sh
PROMOTION_TASK=root-case sbatch --array=1 submit_promotion_battery.sh
```

After those records are healthy, submit the full arms:

```bash
PROMOTION_TASK=bound-parity sbatch --array=1-29%29 submit_promotion_battery.sh
PROMOTION_TASK=demand-point sbatch --array=1-17%17 submit_promotion_battery.sh
PROMOTION_TASK=root-case sbatch --array=1-6%6 submit_promotion_battery.sh
PROMOTION_TASK=strict-repeat sbatch --array=1-2%2 submit_promotion_battery.sh
PROMOTION_TASK=diagnostics sbatch --array=1 submit_promotion_battery.sh
PROMOTION_TASK=calibration sbatch --array=1-2%2 submit_promotion_battery.sh
```

Collect after every task finishes:

```bash
cd code/model
PYTHONPATH=$PWD python -m intergen_housing_fertility_optimized.promotion_collect \
  --run-root ../../output/model/intergen_optimized_promotion_20260719
```

The collector fails nonzero unless all six gates pass. It writes compact CSV
sidecars, both selected-case 15-row target-fit tables, `summary.json`, and
`REPORT.md`.
