# Fixed-coordinate parenthood housing floor diagnostic

Submitted October 1: replay gate **19000409**, followed by array
**19000410_0–3** with an `afterok` dependency. The first independent scheduler
check showed the gate pending on priority and the array pending on dependency.
No numerical result is claimed. The lead independently passed the local mock
and verified the staged runner, incumbent and manifest SHA-256 values exactly.
The unchanged existing calibration and transition jobs were not touched.

This packet tests the physical parenthood housing floor $h_P$ at 2.3, 2.4,
2.5, and 2.6 rooms. It is a fixed-coordinate sensitivity, **not a
recalibration**. The nine other coordinates stay at the normalized v2 best
postchecked point, chain 2 case `0064_nm` ($h_P=2.2992335366442824$, base
loss 29.969678753593804). The selected price 0.7136094701329704 is only a
root starting value; the existing fertility renewal price root remains active.

The sole economic change relative to this incumbent is the proposed physical
floor $h_P$. The old search bound [0.1, 2.3] is extended locally to [0.1,
2.6] so the four evaluations are admissible; the bound is metadata and does
not change utility at a fixed point. Population remains $N_0=1$ and $H_0$ is
derived internally at the actual demand and accepted price. Earnings, entry,
transfers, preferences apart from $h_P$, mortgage and credit, target values,
weights, 120-by-9 grid, and reporting gates retain the normalized v2 contract.
No search, target change, or parameter adoption is implied. The ten scored
moments for ten previously free coordinates do not establish local rank here.

`incumbent.json` records the terminal Torch postcheck source and SHA-256;
`manifest.json` pins the v2 driver, normalization adapter, plan, source pins,
target fingerprint `db60605e...ba1`, and weight fingerprint `2391cd2d...0b0`.
The runner calls the v2 native evaluator directly. The gate first evaluates
the **exact incumbent** under old and extended bounds and requires identical
target tables, all 31 estimates, price, derived $H_0$, and exact native
ROOT/REPEAT reports and 17 plots. Only after that gate passes can the four
points run. Each point writes the full 14-row target fit, 31-row parameter
table, and native 17-plot ROOT/REPEAT diagnostics. `latest_completed.json` and
`best_so_far.json` are written for each completed point. The Slurm launcher
writes a terminal receipt on every exit.

The prior v2 selected postchecks took 169–243 seconds, median 213 seconds.
The six planned native evaluations (two gate replays plus four parallel points)
should take roughly 7–9 minutes for the gate and 3–5 minutes for the array
after queue/startup; both Slurm stages have a hard 30-minute walltime per
worker, one CPU, 24 GiB, and one math thread. There is no automatic retry.

Preparation and launch:

```sh
PYTHONDONTWRITEBYTECODE=1 code/model/.venv/bin/python output/model/fixed_reference_economics_20260928/normalized_floor_extension_v1/run.py --mode mock --out output/model/fixed_reference_economics_20260928/normalized_floor_extension_v1/mock_check
bash output/model/fixed_reference_economics_20260928/normalized_floor_extension_v1/stage_torch.sh
bash output/model/fixed_reference_economics_20260928/normalized_floor_extension_v1/submit_torch.sh
```

The stage script does not submit. The submit script sends one gate job and an
array of four point jobs with `afterok` on the gate. Existing v2 jobs and
results remain untouched. The outcome should be reported as target and model
moments, full loss contributions, and all fixed/free parameter positions,
without interpreting this diagnostic as a fitted calibration.
