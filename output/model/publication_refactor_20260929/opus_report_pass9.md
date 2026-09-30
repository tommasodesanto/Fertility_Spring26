# Opus pass-9 report — native full-GE benchmark (verification only)

Identity: Claude **Opus 5.5** (`claude-opus-5-5`), session `861c5a63-40cb-4bde-ab69-9caec23d0576`. The pass-8 report is kept as `opus_report_pass8.md`. No model, lifecycle or GE execution, jobs or commits. No source rewrites, bypasses or monkey-patches; `run.py` is unchanged and independent of the original sources.

## Implemented: `refactor_lab/native_ge_benchmark.py`

- **Command:** `run --engine original|lab`. It requires `--credit reference`, `--initial-price-factor 1.05`, a fresh empty `NUMBA_CACHE_DIR`, and all thread variables set to 1.
- **Inputs:** both engines use the same bundle and external SHA, the same grid, and the same complete initial public parameters (written as `parameters_initial.json`).
- **Original engine.**
  - Before import, the 13 pinned sources are authenticated against `engine_receipt_pass2_reference.json`, and the frozen `source_manifest` hash against the reference manifest.
  - It then imports the genuine `intergen_eqscale_seq_optimized.solver` and `tools/e5f_stationary_paygo.py`. The frozen calibration or calendar runtime is never imported.
  - After import and after the solve, every imported file under `ROOT/code/model` must match a pass-2 pin or the frozen manifest. Any unpinned or changed file stops the run.
- **Lab engine.** Imported from `PYTHONPATH` (scalar stage or the indexed candidate). Every imported `refactor_lab` file hash is recorded; any original package module in the process is a failure.
- **Solve.** Identical `solve_balanced_initial_equilibrium(..., payroll_tax=P.tau_pay, 1e-9, 1e-6)` for both. Same `CallCounter` (18 at-price calls, 900 s solve stage, SIGALRM backstop); budget stop exits 4.
- **Native gates, required:**
  - strict market convergence;
  - the returned pension `marginal_gate`/`fiscal_gate`;
  - parameter identity (only `native_inherited_distribution_evidence_dir`, `eq_iter` and values derived by the engine's own `bind_initial_balanced_pension` may differ);
  - finite price and distribution.
  - Renewal gap is classified against the reference `adult_entry_gate` with a 1e-6 diagnostic threshold; no psi change.
- **Saved outputs:**
  - every solution and shared array, plus `benchmark.b_grid`, `benchmark.price` and `benchmark.start_price`;
  - the stage pickle, and `parameters_final.json`;
  - a receipt with runtime identity, threads, cache, `PYTHONPATH`, executable and source identities, profiled counts (including payload upgrade vs replay), and timing (inputs, import, profiled solve stage, serialization, in-process total). External wall time and sampled RSS come from `budget_run.py`.
- **`compare-params A B`:** full effective-parameter comparison, excluding only run fields. Arrays use the existing `compare.py`: exact, with nothing missing.
- **Label:** `native_GE_benchmark_only`. It makes no claim about the 113 paths, 14/31 tables or 17 standard plots; those stay on Torch (job 18849543 and later). No supplemental plots are generated.

## Checks actually run (import/authentication/CLI only)

- **Authentication:** all 13 pins verify, and the frozen `source_manifest` hash matches.
- **Import audit:** after importing the genuine original solver and PAYGO, 11 imported ROOT files are authenticated (10 pass-2 pins; the package `__init__.py` via the frozen manifest).
- **Lab audit:** refuses a process that has imported original modules.
- **CLI guards:** refuses factor 1.10 and non-single-thread environments.

## Commands

See the "Native full-GE benchmark" section of the README: `original`, `lab`, and optionally `indexed` (`PYTHONPATH=<indexed_src_pass8>`), each under `budget_run.py` (1200 s, 12 GiB sampled) with a fresh cache. Then:

- `compare.py` (lab or indexed vs original);
- `native_ge_benchmark compare-params`.

## Caveats

- Modules imported lazily during the solve (e.g. `fertility_nested`, `two_shock_choice`) are audited only after the solve. A changed file would be detected after it had already executed; the run then fails, but no earlier stop is possible without an import hook.
- Timing includes identical profiler overhead in both engines.
- ARM results are not compared with the Torch checkpoint here.
