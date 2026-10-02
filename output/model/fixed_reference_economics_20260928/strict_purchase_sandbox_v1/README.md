# Strict purchase origination sandbox

This is a **single fixed-coordinate experiment**, not an adopted model or a calibration search. The reference is normalized v2, chain 2, case `0064_nm`, with all ten parameters fixed. The only economic change is that wealth held at the start of the period, plus proceeds from sale of an old home, must cover the required down payment. Current income cannot qualify a new purchase. The purchase check is \(b+S\geq(1-\phi)Q\), with post-transaction balance \(b+S-Q\geq-\phi Q\). Income still enters the unchanged budget \(b'=R(b+S-Q)+y-c-K\). Interest timing, final debt bound, stayer rules, targets, weights, and normalization remain unchanged.

The two engine copies live at `code/model/experiments/strict_purchase_sandbox/source/`. Each was copied from its original refactor or indexed-credit source, and only the `purchase_income` block in `engine/household.py` changed. `native_purchase_income=True` remains set to retain the existing transaction support map and kernel interface. `manifest.json` pins the original and sandbox household/kernel files, incumbent, target and weight fingerprints. `verify_fixture.py` checks every copied source byte except the intended block and executes a minimal origination/budget fixture without a model solve.

Run `python3 verify_fixture.py` and `bash -n stage_torch.sh launch_torch.sh` before staging. `stage_torch.sh` prepares a new Torch directory without submitting. The staged `launch_torch.sh` requests one CPU, 24 GiB and 20 minutes for one native normalized GE root and repeat; `run.py` requires 14 target rows, 31 parameter rows, 17 standard plots and repeated-report agreement. The launcher uses the frozen runtime and isolated source mounts; it has no search or automatic retry.

## Torch run

Job **19004856** (`strict_purchase_one`) was confirmed **RUNNING** on `cs655`; Slurm reports start **2026-10-01 23:44:58 New York** and a **20-minute** limit (deadline 2026-10-02 00:04:58). It has **one CPU and 24 GiB**. The submission command from the staged directory was:

```bash
sbatch --job-name=strict_purchase_one launch_torch.sh
```

The staged directory is `/scratch/td2248/projects/strict_purchase_sandbox_v1/`. The launcher refuses an existing `results/onecase` directory, so this command documents this submission rather than an in-place retry. There is no completed result or economic comparison yet.

The fixed reference is normalized v2, chain 2, case `0064_nm`. Its ten coordinates, all held fixed during the experiment, are:

| Parameter | Value |
| --- | ---: |
| `beta_annual` | 0.9670843817494936 |
| `chi` | 1.093429912990833 |
| `child_benefit_curvature` | 0.09889818563569602 |
| `first_birth_fixed_cost` | 0.38540954564501306 |
| `h_P` | 2.2992335366442824 |
| `kappa_fert` | 0.126367646227124 |
| `kappa_fert_continuation` | 0.3511943513636155 |
| `psi_child` | 0.17141586769606804 |
| `tenure_choice_kappa` | 0.013862580239215504 |
| `theta0` | 0.10416394115934847 |

The closure keeps initial physical population \(N_0=1\), normalizes completed fertility to 2.1 with the fertility-price root, and internally derives the housing-supply coefficient \(H_0\) for each equilibrium. The incumbent's \(H_0=6.866751098158112\) is a reference value, not a fixed sandbox input. The incumbent price 0.7136094701329704 is only the GE starting price. The financed share remains \(\phi=0.8\); the real interest rate, pension PAYGO rule, targets, weights, entry distribution, grids and all other economic inputs follow the normalized-v2 contract.

The lead verified the staged source files against these local SHA-256 pins:

| Source | Original | Sandbox |
| --- | --- | --- |
| `refactor_lab/engine/household.py` | `2a34f5f28c0c63ca3d24d9e759cd1b5b92aadaddb5fe04795e65abbc0b448082` | `3c67efbb2738dfc1f7a2157e92a00330dfab0ed4ca6d55608ddc90fc8fb829d2` |
| `small_credit_lab/engine/household.py` | `2a34f5f28c0c63ca3d24d9e759cd1b5b92aadaddb5fe04795e65abbc0b448082` | `3c67efbb2738dfc1f7a2157e92a00330dfab0ed4ca6d55608ddc90fc8fb829d2` |
| Both `engine/kernels.py` copies | `fe7d43afdac234af71edd09e1260666e21c309dda83de4adc6a3abff35c3d5d4` | Same bytes |

The target fingerprint is `db60605ef444b747b2ed7f9482b6b8c90ef5e8b1d75c49a3ac1ff6f4275c9ba1`; the weight fingerprint is `2391cd2d4a39a6669a405be34ff14116ad314354353173cb619f7fa7c66043b0`.
