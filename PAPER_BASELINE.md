# September 14 paper baseline

This branch is active `main` development for the model shown in the September
14 presentation. The immutable Monday reference is tag
`paper-baseline-2026-09-14`; 365 of 367 manifest sources match that tag
byte-for-byte. Two main-only builders,
`build_simplified_olg_theory_slides.py` and
`build_e5f_bounded_refinement_plan.py`, are intentionally preserved.

The code is commit `70abd4a8` with the exact corrected `solver.py` used by
job 17655042: SHA-256
`2992412586b81cef3a3e58d92191bb51f54d3f9cc600d7675bbadaed7d1682da`.
The 365 restored source files match the frozen tag byte for byte. The
pension helpers and exhaustive saving solver are included. The snapshot's
experimental choice modules are retained for source identity; the retained
calculation uses sequential household choices.

The reference retains balanced PAYGO pensions, payroll tax 0.179, the original
household birth queue with the 2.1 replacement normalization, no immigration,
and the baseline 1% annual property tax with equal rebates. Preservation of
these choices does not resolve the outstanding economic questions in the ledger.

## Check the reference

From this worktree:

```bash
python3 code/model/tools/check_paper_baseline.py --revision HEAD
```

Use `--revision HEAD` to check the committed source because ongoing local edits
may be present in the `main` checkout. This checks the frozen source and
evidence hashes without solving the model.
The existing fiscal-accounting and exhaustive-saving tests also passed (36
tests). A fresh end-to-end numerical replay has **not** been performed during
this freeze.

## Original recipes and results

`output/model/paper_baseline_sep14/manifest.json` pins the source, recipes and
retained outputs. `initial_recipe/` contains the original initial-state launch
scripts. `native_files.json` maps original cluster paths to archived local
copies under `native/`, including the operative run contracts, scoring code,
and announced-transition driver and dependencies.

The operative initial run contract names `corrected_initial_source_v2`.
The old `preparation.json` and fetched `MANIFEST.md` retain earlier source names
and hashes; they are preserved as historical evidence, not used to select the
solver. Do not execute the archived `run.sh` unchanged: it contains absolute
cluster paths and the original output destination. A replay must use a fresh
output directory and explicitly verify/relocate its dependencies.

`retained_results/` includes the full initial target and parameter tables, the
saved presentation PDF, and the plotted transition data. The saved presentation
transition is the **four-announced-shock** curve from
`announced_original_queue_20260913c`; its first fertility point is
1.8279799986983412. It is a provisional, unconverged transition. Some historical
slide captions and folder names incorrectly describe it as a single permanent
shock. Preserve the actual curve as the reference; do not use those labels to
choose a different run. The later successive-surprise recovery is not this
presentation reference.

## Next work

1. Replay the original initial-state recipe and capture its fully resolved
   parameters and inputs; reconcile the sandbox against that recipe.
2. Establish numerical regression tests against this frozen reference before
   testing new economic specifications.
3. Measure runtime, refactor, and require unchanged outputs and numerical gates
   for changes presented as efficiency improvements.

Keep the reference commit/tag immutable. Active `main` development retains
ongoing author, sandbox and refactoring work; the immutable tag remains the
paper reference.
