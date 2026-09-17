# Fetched files — corrected_initial_template_v6 driver (read-only cluster fetch)

Fetched 2026-09-16 via `scp` from Torch (`ssh torch`), no cluster writes, no jobs touched.

Base cluster path: `/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/final_night_20260913/`

`source_root` for the pension functions (named in `preparation.json`):
`/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/final_night_20260913/corrected_initial_source`
(a 641-file snapshot of the repo at commit `70abd4a8` with one corrected file,
`code/model/intergen_eqscale_seq_optimized/solver.py`; only that one Python file plus the
pension-function file below were fetched, not the full snapshot.)

| Local file | Size (bytes) | SHA-256 | Cluster path |
|---|---|---|---|
| run.sh | 1370 | 666df5f0694d88e67bceb4bf5ca4e764b02e983634ba872b706cd450905f5e37 | batches/final_night_20260913/corrected_initial_template_v6/run.sh |
| run_capped_beta.py | 20587 | 8b5ed3c50ee32805fea8147d6f2f617b66ad239218e8a90c3019e31593e26ee3 | batches/final_night_20260913/corrected_initial_template_v6/run_capped_beta.py |
| run_profile.py | 29859 | a22f2f9044f1d7ea737d45c7ee57b94c263ba302f9bce4c7cff4b6fb35d2df1f | batches/final_night_20260913/corrected_initial_template_v6/run_profile.py |
| preparation.json | 1807 | 3a60a7e9b551243998608519f65f3726773e611c719798c9d482796da06543a4 | batches/final_night_20260913/corrected_initial_template_v6/preparation.json |
| proposal.json | 454 | a92f642a9b71e7a29de79a3e1cd604163825b8c68efd8ab95e10d43f02db68c7 | batches/final_night_20260913/corrected_initial_template_v6/proposal.json |
| submission_v6.json | 1812 | 8ac2c92a683edd30ebe06b71f4f6c147185dc356c6a409c659c99c2603a2b82f | batches/final_night_20260913/corrected_initial_template_v6/submission_v6.json |
| source_objective_version.json | 837 | 6722c900dacfd6498b5efab0fcfdd391e53ef1b0e7bca8f126811526707eed48 | batches/final_night_20260913/corrected_initial_template_v6/source_objective_version.json |
| run_e5f_joint_rebated_initial_scored.py | 3420 | 9da35b55466d74dc10a6a85ff91dabe887a807aab68417eabf63a9d7a42ca967 | batches/final_night_20260913/joint_source/run_e5f_joint_rebated_initial_scored.py |
| run_e5f_rebated_initial_overnight.py | 19704 | d9aa97b890442d45971ec622b4b41687da10ffa26eedaf0198f352c7e6ecb790 | batches/final_night_20260913/initial_v2/run_e5f_rebated_initial_overnight.py |
| run_e5f_joint_rebated_initial_probe.py | 12798 | 9eee3bca39f2a98f4a58cf18196d695b6a9db3e93ef18f8eaa2bf4dbb1243bbb | batches/final_night_20260913/joint_source/run_e5f_joint_rebated_initial_probe.py |
| e5f_stationary_paygo.py | 7664 | d723bae3e2a6d912f1095a5aa94a1780eb9ad6a43babdd45659bea31c680a6f9 | batches/final_night_20260913/corrected_initial_source/code/model/tools/e5f_stationary_paygo.py |

Notes:
- `run.sh` invokes `run_e5f_joint_rebated_initial_scored.py` with `--helper`
  (`run_e5f_rebated_initial_overnight.py`), `--joint`
  (`run_e5f_joint_rebated_initial_probe.py`), and `--template` pointing at this folder.
  The `--helper-sha256` / `--joint-sha256` args in `run.sh` match the fetched files' hashes
  above, confirming these are the exact source files the run used.
- `joint_initial_solution` is defined in `run_e5f_joint_rebated_initial_probe.py` (line 106).
- `bind_initial_balanced_pension` (line 51) and `certify_initial_pension` (line 68) are
  defined in `e5f_stationary_paygo.py` under the `corrected_initial_source` snapshot, not
  inside the template folder itself — fetched from `source_root` per the task instructions,
  rather than pulling the full 641-file snapshot.
- `run_capped_beta.py` is the template's own controller (imports `run_profile.py` as
  `support`); it orchestrates frozen scored-candidate evaluations and does not itself define
  the pension functions.
- Skipped: `terminal.pkl.gz` / other checkpoints and any file over 20 MB (none in this
  template folder exceeded that; the `__pycache__/*.pyc` files were also skipped as
  non-source build artifacts).
