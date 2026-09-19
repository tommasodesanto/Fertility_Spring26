# Local laptop setup

Verified on September 19, 2026, macOS 26.6 / Apple Silicon. This is runtime
readiness evidence, not certification of a full calibration or policy transition.

## Python model

Use `code/model/.venv/bin/python` from the repository root. System `python3`
is 3.9.6 and is too old. The migrated environment depended on an Intel-only
Miniconda executable and failed with `bad CPU type in executable`.

The replacement is native CPython 3.10.21, managed by `~/.local/bin/uv`.
NumPy 1.24.3, Numba 0.58.1, llvmlite 0.41.1, pandas 1.5.3 and matplotlib
3.7.1 retain the versions found in the migrated installation. Other installed
dependencies are pinned in `code/model/requirements-laptop.txt`; these pins
record this tested laptop environment, not a recovered cluster lockfile.

To recreate a missing environment (preserve an existing one first):

```sh
~/.local/bin/uv python install 3.10.21
~/.local/bin/uv venv --python 3.10.21 code/model/.venv
~/.local/bin/uv pip install --python code/model/.venv/bin/python -r code/model/requirements-laptop.txt
```

The old environment is preserved outside Git at
`~/Library/Application Support/FertilityMigration/20260919/model-venv-intel`.
Its original Intel Miniconda installation was not changed.

Run the focused regression checks from the repository root:

```sh
PYTHONPATH=code/model:code/model/tools NUMBA_DISABLE_JIT=1 OPENBLAS_NUM_THREADS=1 \
  code/model/.venv/bin/python -m pytest -q \
  code/model/tools/test_calendar_policy_reuse.py \
  code/model/tools/test_e5f_matched_pf_smoke.py \
  code/model/intergen_eqscale_seq_optimized/tests/test_calendar_time_continuation.py

PYTHONPATH=code/model:code/model/tools NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  code/model/.venv/bin/python -m pytest -q \
  code/model/intergen_eqscale_seq_optimized/tests/test_calendar_time_continuation.py \
  code/model/intergen_eqscale_seq_optimized/tests/test_optimized_port.py \
  code/model/tools/test_e5f_pf_choice_ownership.py \
  code/model/tools/test_e5f_simple_fertility_nest.py
```

Results: 19 pure checks and 29 checks with JIT enabled passed (the four calendar
checks appear in both groups). These exercise tiny-grid model solutions,
calendar continuation, price-search consistency, policy reuse and choice
probability accounting. A stale mock in the fertility-nest test was updated to
accept and verify the current kernel's stayer-value argument; production code
was unchanged. No full-size calibration, historical fit or policy run was made.

Native Numba nopython compilation, PDF plotting, and pandas Stata/Parquet
write/read round trips also passed. Logs, runtime versions and smoke artifacts
are under `output/setup/laptop_20260919/` (local generated files).

## Other workflows

| Workflow | Verified status |
| --- | --- |
| Stata | Native StataMP 17 batch execution passed a generated-data assertion. Binary: `/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp`. Not on the shell PATH. |
| LaTeX | `latexmk -pdf` compiled a temporary copy of `latex/JMP_DS_draft` to four pages. Source and author PDF untouched. TeX Live 2021; an empty-bibliography warning remains. |
| R | Blocked: `/usr/local/bin/R` and `Rscript` resolve to Intel R 4.5 and fail on this machine. Install native Apple Silicon R, then restore/check the packages required by each empirical driver. |
| Torch | SSH configuration exists and the login host is reachable, but authentication is rejected. Renew interactive access with `ssh torch`, then verify `bash code/cluster/torch.sh status`. No job submitted. |

The ATTOM assessor `.dta`, AHS raw data, PSID outputs, MMS family-size outputs,
and mortgage-policy birth data were found locally. Follow-up verification
corrected the initial inventory: CPS `code/data/cps_fertility/cache/jun24pub.csv`
is present (145,653,216 bytes), and all 37 NCHS natality files listed in
`code/data/nchs_natality_timing/first_birth_counts_manifest.csv` exist at their
external archive paths with matching byte sizes. Their hashes were not rerun.
Only the original ATTOM shard directory referenced by `merge_to_stata.py`,
`/Users/tommasodesanto/Desktop/019e6b85-e6d0-79f3-ac21-a5dc416c8dfa`, was confirmed
absent. The merged ATTOM data remain available; the missing source folder
affects rebuilding that merge. This was a targeted availability check, not a
full data inventory or empirical reproduction.

The checkout had substantial pre-existing edits. They were preserved. The
September 19 generated memory snapshot incorrectly described it as clean;
live Git status was used for this setup work.

## Full-grid timing replay, September 19

One saved-price stationary replay of the September 14 paper baseline completed
on the Apple M5 Pro (48 GiB RAM) in **16.52084 seconds**. It recomputed the full
household lifecycle and stationary distribution at the saved equilibrium price:
17 ages, 120 wealth nodes, all original income/family/tenure states, one thread.
It did not search for prices, recalibrate parameters or solve a transition.

The matching frozen source is `tmp/paper_baseline_sep14`, checked against the
paper manifest; the current-main solver has a different hash and was not used
for this historical timing comparison. The driver reads the September 17
Torch checkpoint using compatibility mappings for NumPy/Python pickle namespace
renames; numerical data and parameters are unchanged.

| Phase | M5 Pro, new measurement | Torch, saved September 17 measurement |
| --- | ---: | ---: |
| Household backward solution | 11.19233 s | 13.51177 s |
| Stationary distribution and statistics | 5.32488 s | 7.79160 s |
| Sum of those phases | 16.51721 s | 21.30337 s |

The laptop used an initially empty Numba cache; its timing includes compilation.
The Torch phase times come from the final solve inside a warmed 17-evaluation
initial root, not a fresh benchmark job. Runtime versions differ, and the
cluster CPU model has not been recovered. The observed laptop phase total is
about 22.5% lower, but this is not a controlled hardware speedup estimate.
The 384-second Torch root time covers multiple household solves and must not
be compared with this one-solve laptop timing. No matched old-Mac measurement
is available yet; old-Mac access and renewed Torch authentication are pending.

All ten policy/value/distribution array comparisons pass at `rtol=atol=1e-9`;
the largest absolute difference is below `1e-12`. Housing relative residual
is `4.5167e-7`, below the unchanged `2.5e-5` gate, and stationary pension checks
pass. The standard 17 diagnostic PNGs were generated; the market plot and
age-30 policy panel were visually inspected. Existing policy kinks are retained,
not diagnosed or changed by this portability benchmark.

Receipts and diagnostics: `output/model/laptop_benchmark_20260919/`.
To run the same single solve on another machine, use its matching source and
checkpoint copies and a new output directory:

```sh
code/model/.venv/bin/python code/model/tools/benchmark_saved_stationary.py \
  --source-root tmp/paper_baseline_sep14 \
  --checkpoint output/model/paper_baseline_sep14/replay_20260917/native_output/raw/repetition_01/initial_state.pkl.gz \
  --manifest output/model/paper_baseline_sep14/manifest.json \
  --output output/model/laptop_benchmark_new_machine
```

The driver pins one thread, writes provenance before solving, refuses an existing
output folder, and stops after 15 minutes. Submit it through the established
Slurm workflow on Torch, not on a login node.
