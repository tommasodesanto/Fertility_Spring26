# Broader physical parenthood floor calibration v1

Prepared experimental full calibration; no model solves, staging or submission by this packet's worker. Reference: normalized v2 chain2/case0064_nm, authenticated incumbent loss 29.969678753593804 (incumbent.json is authoritative). The sole economic search change is h_P upper bound 2.3 to 2.6; lower bound remains 0.1. All ten coordinates remain free. Target/weight fingerprints, other bounds, D=0 nonnegative-mean entrants, 120×9 grid, N0=1, derived H0, renewal-price root and native guards are unchanged. Identification rank remains unverified.

24 deterministic new starts: six per anchor h_P=2.2992335366442824,2.4,2.5,2.6. Each group has one exact anchor and five joint neighbors. Python Random seed 20261003, independent uniform perturbations using the v2 modest halfwidths recorded in plan.json; clipping only at unchanged bounds except extended h_P. Raw and actual displacement and clipping are recorded. Every start uses authenticated incumbent price 0.7136094701329704. The optimizer retains the v2 25% simplex steps.

Each chain: three hours from actual start, 150 objective-call cap, 900 seconds reserved for a separate-process selected-point native postcheck, one CPU/24 GiB/one thread. Budget interruption is provisional, never optimizer convergence. No automatic retry or extension. Search uses existing exploration adapter; selected ROOT/REPEAT native report must pass 14 targets, 31 parameters and 17 exact plot repeats.

Deployment authenticates the successful companion normalized_floor_extension_v2 old/extended-bound gate, its complete manifest, v2 source-pins bytes and inherited source files, byte-identical incumbent and all 14 target rows/31 parameter estimates; only h_P bound and resulting near-bound metadata may differ. Gate Slurm terminal receipt must match explicit FLOOR_GATE_JOB with exit 0. No duplicate native gate is launched.

Local zero-solve check:
`PYTHONDONTWRITEBYTECODE=1 code/model/.venv/bin/python output/model/fixed_reference_economics_20260928/normalized_floor_calibration_v1/test_configuration.py`

After lead review, stage only with deployment/stage_torch.sh. Submission requires completed successful companion gate:
`ssh torch 'FLOOR_GATE_JOB=JOB_ID bash /scratch/td2248/projects/normalized_floor_calibration_v1/submit_torch.sh'`
The submitter checks successful Slurm state and sets afterok; each chain independently checks authenticated gate evidence before any solve. deployment/controller.diff records the controller changes; normalized_objective.py and source_pins.json are verbatim copies of v2. Do not reuse a result directory or bypass the gate.
