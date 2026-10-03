# Legacy model archive — October 3, 2026

This archive preserves older model experiments and results whose implementation is historical. The frozen observer authentication
inventory still requires original access paths for 19 legacy Python source files. It does not replace the canonical stationary-GE package
in `code/model/production/`, the retained transition stack, authenticated oracle
and observer inputs, active experiments, or the frozen September 14 reference.
No numerical run, recalibration, or change in economics was performed.

## Retained contents

| Former location under `code/model/` | Archive location | Role |
|---|---|---|
| `intergen_housing_fertility_howard_test/` | Same basename here | Unpromoted June Howard experiment; 12 files |
| `intergen_surrogate_calibration/` | Same basename here | June surrogate methods prototype; 9 files |
| `benchmarks/` | Same basename here | Historical May benchmark JSONs, figures and README; 62 files |
| `PLAN.md` | `PLAN.md` | Superseded May 7 port plan |
| `run_intergen_model.py` | `run_intergen_model.py` | Historical June one-shot diagnostic implementation |

`manifest.json` records every old/new path and SHA-256 hash. All 85 moved files
were checked: 84 are byte-identical. The runnable historical runner changes only
`model_dir` to the maintained `code/model` location so its mechanics-packet
subprocess still resolves. Its exact pre-move bytes are retained in
`source_snapshot/run_intergen_model.py`; the old `code/model/run_intergen_model.py`
retains its exact original regular-file bytes because frozen source authentication
pins its hash and its original `__file__` depth determines repository paths. Edit
historical archive settings in the archived runnable implementation. The original
runner remains an authenticated compatibility source, rather than a new entrypoint.

The four historical benchmark helpers remain under `code/model/tools/` with
exact original bytes because the frozen authentication inventory pins them too.
Their original `benchmarks/` defaults therefore remain historical pointers.
Use explicit CLI paths: `grad_descent_report.py --gd-json` should select
`calibration_archive/model_legacy_20261003/benchmarks/grad_descent_bench.json`;
`gradient_descent.py --out`, `mini_parameter_solve.py --json`, and
`plot_lifecycle.py --out` should select a new path under
`output/model/legacy_benchmarks/`. Their interim default-path edits were reverted
to preserve frozen pins. Exact helper snapshots remain retained here.

## Historical import and command context

From the repository root, include both the archive and maintained model sources:

```sh
export PYTHONPATH="$PWD/calibration_archive/model_legacy_20261003:$PWD/code/model${PYTHONPATH:+:$PYTHONPATH}"
```

The surrogate prototype still imports the retained June package
`intergen_housing_fertility` and retains its original absolute `REPO_ROOT` in
`data.py`. That historical machine-specific setting must be reviewed when
running elsewhere. No training data or outputs were moved. Howard command
default output paths are relative to the working directory; use the original
`code/model` working directory or pass an explicit `--outdir` when inspecting
historical commands. Do not run commands merely to validate this move.

The original prototype READMEs and plan are historical evidence: their old
commands and claims are preserved, not promoted to current documentation.

## Verification and limits

Verification used SHA-256 checks for all moved files, syntax-only AST parsing,
imports of the two archived packages and selected helper modules, an imported
compatibility runner, and a mocked subprocess interception that checked its
mechanics-packet path, source path and working directory. Interim helper parser defaults
were evaluated without importing their solvers, then reverted because frozen
authentication pins the original helper bytes. Native solve count: **zero**.
This proves source preservation and the checked path/import behavior; it does
not certify historical numerics or replay results.

## Authenticated compatibility paths

The initial import-only review missed source authentication references in
`output/model/overnight_calibration_20260928/contract_v1/source_manifest.json`.
Production reporting authenticates that frozen inventory through
`output/model/fertility_identification_20260928/contract_v1/contract.json` and
`fixed_reference_manifest.json`. It pins ten Howard Python files, eight surrogate
Python files, and the original runner: 19 legacy source files in total.

Relative symlinks at `code/model/intergen_housing_fertility_howard_test` and
`code/model/intergen_surrogate_calibration` restore access to archived source.
All 12 Howard and nine surrogate archive files remain byte-identical; all 19
frozen Python pins match through the original paths. The original regular runner
was restored from its exact snapshot, preserving both its hash and path depth.
The archive therefore holds the historical implementations while the two old
package paths and original runner remain authenticated compatibility dependencies.
No frozen manifests, hashes, bundles, or numerical code were changed. See
`authenticated_compatibility` in `manifest.json` for exact identity and receipts.

The final cleanup-owned inventory comparison found 24 frozen Python pins: the
19 legacy package/runner files, `audit_intergen_housing_ladder.py`, and all four
benchmark helpers. All five helper files were restored from `590c9f8b^` after
checking those bytes against the frozen digests. All 24 cleanup-owned pinned
paths now match exactly. Historical explanatory pointers and default paths in
authenticated helpers remain unchanged; this README gives archive navigation
and explicit benchmark path overrides.
