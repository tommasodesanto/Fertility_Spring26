# E5F refactor baseline inventory

The proposed behavior reference is the exact last certified recovery runtime recorded on 2026-09-15, described by `output/model/e5f_original_queue_20260913a/long_successive_refit/recovery_source_manifest_20260915.json`. Its status is **outstanding**: the recovery source directory and all scratch inputs were not found in the inspected local paths/history, and the cluster is intentionally inaccessible for this inventory.

The active chain is `code/cluster/run_e5f_long_successive_refit.py`, prepared by `code/cluster/prepare_e5f_long_successive_refit.py`, using the copied scaled-step root and Toeplitz Jacobian helpers under `code/model/tools/`. The recorded and current SHA-256 hashes match for the driver (`d53ddbdc…`), preparer (`1486bed9…`), and Toeplitz helper (`94c44836…`). The scaled-step helper does not match: certified recovery recorded `ffcb2334…`, whereas the current local file is `c26d312a…`. Git history identifies the current family as the later diagnostic change committed in `25a8b6bd` (“Trimmed-score acceptance…”), after `da875408` (“loosened-gate mode”). The exact `ffcb…` bytes are unavailable in local scratch and cannot be reconstructed from the tracked history inspected here. Therefore the current helper is not silently accepted as the baseline.

Locally available supporting model files include `code/model/tools/e5f_original_queue_terminal.py`, `code/model/tools/run_e5f_transition_calibration.py`, and `code/model/tools/run_e5f_successive_surprises_overnight.py`; `e5f_social_security.py` is absent. NumPy and Matplotlib are importable. The recovery manifest’s scratch `spec.json`, empirical blocks, source manifest, root receipts, rows, terminal pickles, and batch scripts were not found at their `/scratch/td2248/...` paths. The local retained artifacts and manifests remain available under `output/model/e5f_original_queue_20260913a/`.

Reproducible inventory commands are:

```bash
git status -sb
sha256sum code/cluster/run_e5f_long_successive_refit.py code/cluster/prepare_e5f_long_successive_refit.py code/model/tools/e5f_ssj_scaled_step_root.py code/model/tools/e5f_ssj_toeplitz_jacobian.py
git log --all --oneline -- code/model/tools/e5f_ssj_scaled_step_root.py
```

The complete pre-existing working-tree status is saved in `baseline/git_status_full.txt`. No model solve, test, installation, cluster call, source edit, or status-file edit was made.
