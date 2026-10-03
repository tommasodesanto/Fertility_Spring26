# Alternative wealth target: isolated Torch search

The isolated stage passed host/container source checks and a chain-9 zero-solve
preflight. Torch smoke array **19111315**, task 9 only, was submitted with the
90-minute budget; production remains unsubmitted pending native smoke review.
The stage and launch receipts are in the new run packet's `deployment/` folder.

This package stages the passed October 2 soft-timing source archive and three
explicit new inputs: `cluster_calibrate.py`, the ten-start plan, and the verified
new-wealth chain-02 `completed.json`. It never uploads the dirty working tree.
The only new search control is the annual discount-factor bound
`[0.930, 0.990]`; the post-interest timing, target 4.45838713455674, numerical
weight 7.595098472533724, and remaining economics match the local new-wealth
experiment. This search is experimental and does not adopt a paper calibration.

The isolated Torch root is
`/scratch/td2248/projects/alternative_wealth_calibration_20261003_v1`.
The commands that staged this exact run from the repository root were:

```bash
bash code/cluster/alternative_wealth_calibration/stage_torch.sh
ssh torch 'cd /scratch/td2248/projects/alternative_wealth_calibration_20261003_v1 && sbatch --parsable --array=9 --time=01:30:00 --export=ALL,WEALTH_RUN_MODE=smoke launch_torch.sh'
```

Staging checks every source hash on the host and in Apptainer and initializes
chain 9, whose starting discount factor is 0.930, without a model solve. The
smoke runs two actual optimizer calls and a fresh-child full native check on one
core with 24 GiB and a 90-minute wall time including a 1,800-second native
reserve. Save the returned smoke job ID. Review `results/smoke_alternative_chain_9`
and run `verify_smoke_gate.py` before production. Production is a separate,
lead-reviewed `submit_torch.sh` action; it submits exactly array `0-9%10`,
one core and 24 GiB per chain, six hours per chain including the 1,800-second
native reserve, with at most 500 objective calls per chain. No fallback,
automatic retry, or other cohorts are configured.

Monitor only the new job with its returned ID:

```bash
ssh torch 'squeue -j JOB_ID -o "%.18i %.2t %.10M %.10l %.30j"'
ssh torch 'cd /scratch/td2248/projects/alternative_wealth_calibration_20261003_v1 && find results -name heartbeat.json -maxdepth 4 -print'
```

Each chain writes `heartbeat.json`, `latest_completed.json`, and
`best_so_far.json` after completed cases. `collect_torch.py --mode smoke` or
`--mode production` with `--root` and `--fetch` copies and validates terminal
reports, requiring the 14-row new-target fit, 31 parameter records, and 17
standard diagnostic plots. A numerically verified result remains provisional
if its optimizer has not converged.
