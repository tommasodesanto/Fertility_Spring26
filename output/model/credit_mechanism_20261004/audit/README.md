# Saved native budget and purchase audit

The native checks remain **uncomputed** for both saved cases. A single fresh,
absolute-path initialization for `phi_095_run1` stopped before policy reconstruction
or evaluation, after 1.830 seconds with 1.346 GiB peak resident memory. No Bellman
or stationary solve occurred. The `price_110_run1` check was not repeated because
it uses the same blocked authenticated factory and source contract.

The exact missing file is:

`/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/calibration_archive/model_legacy_20261003/intergen_housing_fertility_howard_test/__init__.py`

The initial working-tree status lists this tracked archive file as deleted.
`production.reporting.build_context` reaches `authenticate_frozen`, then
`e5f_evening_calibration_runtime.setup -> verify_sources`, which resolves the
historical compatibility path into that archive and fails while hashing it.
Absolute output paths therefore do not resolve the first initialization failure.
The previous second failure was the separate no-overwrite guard on an already
created `cached_audits/runtime_auth` directory.

`summary.json` and `audit_status.csv` pin both saved array and executed-parameter
files against each completed-case receipt, and verify all 31 source hashes in
each lead-approved reached-runtime contract. Both identities pass. The original
production snapshot is separately pinned by SHA-256 in the actual attempt
receipt. No source was recovered, no authentication was bypassed, and no economic
input, target, entry, fiscal object, manuscript, or gate was changed.

## Evidence and reproduction

- `phi095_absolute_run1/receipt.json`: exact missing filename and full traceback.
- `summary.json`: two-case identities and uncomputed status.
- `audit_status.csv`: compact two-case status table.
- `audit_saved_native.py`: absolute-path cached audit using canonical
  `production.reporting.build_context`, `policy_from_solution`,
  `reconstruct_stationary_pre_fertility`, `evaluate_period`, `dated_budget` and
  `audit_purchase_accounting`. It retains the native mass threshold 2e-10 and
  occupied transaction-wealth threshold 1e-9, checks the solve counter remains
  zero, and refuses an existing output directory. This script has not reached
  the checks after initialization; it is not a completed audit certificate.

The next decision belongs to the lead: recover the exact authenticated historical
source under a separately authorized scope, or leave the budget/purchase checks
uncomputed. Do not rerun into the existing attempt directory.
