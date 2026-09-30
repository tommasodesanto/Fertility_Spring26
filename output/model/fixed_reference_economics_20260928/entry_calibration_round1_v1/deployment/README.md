# Compact deployment

The launcher uses the immutable `/scratch/td2248/projects/grid_resolution_credit053_v2` source stage and frozen September 28 project. It checks all 81 base-source hashes, every new archive-source hash, and both full-grid input file hashes before any native evaluation. The 125 MB full-grid input arrays remain in the existing read-only remote input folder; they are not archived or uploaded. Compact Python, shell, JSON and CSV files named by the final source pins are included, so seed and imported-adapter receipts are available in the frozen container namespace.

Array tasks 0, 1 and 2 correspond to `empirical_credit_120x9`, `nonnegative_mean_120x9` and `nonnegative_mean_160x15`. Each receives one CPU, 24 GiB and exactly four hours from launcher entry, including authentication, its mocked preflight, native repeated baseline, search, final repeats and reporting. Threads are fixed to one. Output folders are created atomically and never reused; failed jobs are not restarted. The EXIT/TERM handler writes `launcher_terminal.json`. An uncatchable SIGKILL cannot execute a shell trap: scheduler accounting is the independent evidence for that case.

After the lead approves final sources and tests, build locally:

```sh
code/model/.venv/bin/python output/model/fixed_reference_economics_20260928/entry_calibration_round1_v1/deployment/build_stage.py
```

Upload/extract only into a freshly created `/scratch/td2248/projects/entry_calibration_round1_v1`; verify the archive against `stage_receipt.json` and every extracted inventory entry. Create its `logs` directory. Do not overwrite an existing remote source or results folder.

Lead-authorized preparatory zero-solve array:

```sh
sbatch --job-name=entry_round1_preflight --time=00:10:00 --export=ALL,CAL_PREFLIGHT_ONLY=1 /scratch/td2248/projects/entry_calibration_round1_v1/launch_torch.sh
```

Each preparatory lane runs the actual test file (300-second cap), then its exact `--mode preflight --lane ...` CLI (300-second cap). It writes under `preflight_validation`, separate from main outputs. Passing this does not certify the native equilibrium; each main lane must pass the native baseline and repeat before search.

After all three preflights and lead gate pass:

```sh
sbatch /scratch/td2248/projects/entry_calibration_round1_v1/launch_torch.sh
```

No jobs have been submitted by this deployment preparation.
