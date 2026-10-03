# Legacy model archive — October 3, 2026

This archive preserves older model experiments and results that have no incoming
active runtime imports. It does not replace the canonical stationary-GE package
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
is a small compatibility forwarder that preserves the historical execution guard
and exported Spyder globals. Edit historical settings in the archived implementation.

The four historical benchmark helpers remain under `code/model/tools/`. Their
numerical code is unchanged. `grad_descent_report.py` reads the archived default
`benchmarks/grad_descent_bench.json`; future default outputs from
`gradient_descent.py`, `mini_parameter_solve.py` and `plot_lifecycle.py` go to
`output/model/legacy_benchmarks/`. Explicit CLI overrides still work. Exact helper
snapshots and path-change receipts are retained here; new outputs must not
replace archived evidence.

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
mechanics-packet path, source path and working directory. Helper parser defaults
were evaluated without importing their solvers. Native solve count: **zero**.
This proves source preservation and the checked path/import behavior; it does
not certify historical numerics or replay results.
