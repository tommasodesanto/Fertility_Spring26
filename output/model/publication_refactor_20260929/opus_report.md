# Opus pass-10 report — reviewed test package

Identity: Claude **Opus 5.5** (`claude-opus-5-5`), session `861c5a63-40cb-4bde-ab69-9caec23d0576`. The pass-9 report is kept as `opus_report_pass9.md`. No model, lifecycle or GE runs, jobs or commits. No original-source or historical receipt edits.

## Diff summary (`code/model/refactor_lab/`)

**A. Promotion — exactly one numerical change**
- `engine/kernels.py` is now the Torch-tested `indexed_src_gridfix` kernel, byte-identical (`fe7d43af…`).
- The unchanged `make_indexed_stage.py` reproduced `kernels.py`, `transform_receipt.json` and `indexed_saving.diff` byte-for-byte from the preserved scalar source before promotion (input kernel SHA `19dceb70…`). Every other engine file is identical to the tested stage.
- Added `engine/transform_receipt.json`, `engine/indexed_saving.diff` and `engine/promotion_receipt.json` (`efc78ea9…`). The promotion receipt:
  - pins the chain: materialized original kernels `19dceb70…` = transform before; promoted = transform after = tested stage; `saving.py` `25e3fa12…`; generator `dd615832…`;
  - pins the evidence: FP 18850069 `completed.json` `3c52051f…`; native pair summary `4359a9c8…` and comparison `291721ea…`; full GE 18851943 pending;
  - status: `accepted_for_fixed_reference_test_package_only`, not an active-production replacement.
- **Tests.** `test_engine_provenance` now checks the explicit chain:
  - before-hash = materialized output, and after-hash = live = tested stage;
  - the new block = `saving.exhaustive_saving_indexed` (renamed) + 4 verbatim helpers;
  - the replaced segment = the original scalar;
  - all other kernels definitions present, and all other transform-listed files unchanged.

  The indexed-kernel test now runs on the live engine by default, against the ORIGINAL scalar oracle (boundary cases + 40 model columns). An external stage remains optional via env vars.

**B. Layout.** Verification-only machinery moved into `verification/` (one folder, no wrapper copies):
- `acceptance_oracle`, `baseline_identity`, `compare`, `budget_run`, `callcount`, `native_ge_benchmark`, `materialize`, `apply_split`, `make_indexed_stage`, `saving`;
- `verify.sh` and `verify_torch.sh`;
- `split_plan.json` and the pass-2 source pins.

Imports and module commands are updated; `ROOT` depth in `verify.sh` is fixed. The root keeps `__init__`, `inputs`, `credit`, `export_inputs`, `run`, `README`, `requirements.txt`, `engine/` and `tests/`. Normal numerical execution imports no legacy engine; the equilibrium driver imports lightweight `verification.callcount`, whose instrumentation is inactive without `--count-calls`. Stale "Torch only" docstrings are fixed.

**C. README.** Rewritten as a user guide:
- scope; one run command with the bundle pin; credit modes, including the D=0 infeasible-entrant warning and the need for an author decision;
- module map and provenance chain; verification commands and the frozen-observer boundary;
- evidence: FP certificate; native pair 152.65 s vs 99.55 s, both 4 household and 5 KFE calls with a payload upgrade, 90 arrays exact; renewal 1.70e-6 against the 7.92e-7 reference (the 1e-6 threshold is a diagnostic difference, not an absolute gate); full GE pending.

## Tests run (local, venv313, one thread)

- **Component suite:** `budget_run --seconds 180 --max-rss-gib 12` around `pytest refactor_lab/tests`, with `REFACTOR_EXPORT=local_export_v1`, bundle SHA `427e67a3…`, the physical checkpoint and the reviewed overlay. **27 passed in 7.05 s**; 8.13 s wall, 1.25 GiB sampled RSS. Evidence is in `pass10_component_tests/`.
- **Driver smoke:** `verification/verify.sh` smoke in local and torch modes passes (`pass10_driver_smoke/`).
- **Dry acceptance:** resolves the tests path and `REFACTOR_ROOT` correctly.
- **CLI:** `--help` works for `run`, `export_inputs` and `verification.{acceptance_oracle, compare, native_ge_benchmark, budget_run}`.
- **Import check:** model imports were checked for `intergen_*` dependencies. `run.py` imports lightweight `verification.callcount`; its instrumentation is inactive unless `--count-calls` is supplied.

## Pending / questions

- Torch full GE/reporting 18851943 is pending; the lead will finalize.
- `transform_receipt.json` keeps its generated status string (`experimental_until_…`). The acceptance decision lives in `promotion_receipt.json`, so tested bytes are not rewritten.
- The equilibrium driver imports `verification.callcount`; instrumentation is optional and off by default. It does not change the numerical model.

## Lead closeout after pass 10

Full Torch GE/reporting18851943 completed successfully. Final receipts are
in `ge_pair_18851943/`; README and promotion metadata now record the result.
Sol's provenance, import-description and launcher findings were corrected
without numerical edits. The source wrapper now uses explicit `LAB_SRC`,
with a relocated-script test. The normal entry remains reference mode by
explicit selection; corrected-zero GE remains blocked by entrant feasibility.
