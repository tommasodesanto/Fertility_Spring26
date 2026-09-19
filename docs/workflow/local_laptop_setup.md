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
| R | Repaired September 19: native R 4.6.1, 33 project packages load, 28 existing NCHS tests and the empirical runtime smoke pass. See the R installation receipt below. |
| Torch | Initial authentication failure was resolved by the author later on September 19. One compute-node benchmark subsequently completed; see below. |

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

## Native R installation receipt

On the author's request, installed the official [CRAN R 4.6.1 Apple Silicon
package](https://cran.r-project.org/bin/macosx/). The downloaded installer's
SHA-1 matched CRAN (`fc9f4ada15589e8e037b9bf05563d21e97181635`), and
`pkgutil --check-signature` confirmed the Apple-trusted Simon Urbanek installer
signature and notarization. macOS administrator authentication authorized the
system installation. The old framework versions 4.2 and 4.5-x86_64 remain;
`R`, `Rscript`, and `/Applications/R.app` now use native ARM R 4.6.1.

Packages were installed as CRAN macOS binaries into the separate user library
`~/Library/R/arm64/4.6/library`, without copying Intel binary packages or changing
the user's `.Renviron`. All 33 packages checked from the empirical/data scripts
load successfully, including `data.table`, `haven`, `fixest`, `ggplot2`, `sf`,
`tidycensus`, `ipumsr`, `modelsummary`, and `tidyverse`. The exact package/version
receipt is `output/setup/laptop_20260919/r/packages.csv`. This is a runtime
snapshot, not a historical package-version lock or full empirical reproduction.

Verification:

- `Rscript code/data/nchs_natality_timing/test_first_birth_timing_targets.R`:
  all 28 existing deterministic checks pass.
- Native architecture assertion; fixed-effects estimate agrees with a dummy-
  variable regression to `1e-10`; clustered standard errors compute successfully.
- Stata write/read numeric round trip, PNG graphics, `sf` coordinate transform,
  and a five-row read of the existing ATTOM Stata file pass. The initial smoke
  compared Stata display-format attributes as well as numeric data; the corrected
  numeric comparison passes with tolerance `1e-12`.

Installation script, smoke script, logs and session details are under
`output/setup/laptop_20260919/r/`. No production data builder was rerun and no
empirical output was overwritten. MATLAB installation remains author-owned.

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
is available yet; old-Mac access is pending. Torch access was subsequently
restored and a fresh matched run completed, as recorded below.

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

### Fresh Torch comparison after login renewal

Job `18029672` completed successfully on `cs602`, an Intel Xeon Platinum 8592+
node, with one allocated CPU and an empty Numba cache. Source-package hashes
and the input checkpoint SHA-256 match the laptop receipt exactly. Every array
comparison and the housing/pension gates passed.

| Phase | M5 Pro | Fresh Torch run |
| --- | ---: | ---: |
| Household backward solution | 11.19233 s | 22.66038 s |
| Stationary distribution and statistics | 5.32488 s | 7.63218 s |
| Full solve call | 16.52084 s | 30.29908 s |

The observed solve-call ratio is 1.834 (45.5% less elapsed time on the laptop).
Both are single cold-cache observations, not repeated warm-speed estimates.
The installed runtimes differ: laptop Python 3.10.21 / NumPy 1.24.3 / Numba
0.58.1; Torch Python 3.13.5 / NumPy 2.1.3 / Numba 0.61.0. This compares the
working installations, not CPU hardware alone. Diagnostics and checkpoint
loading are outside the solve timer on both machines.

Torch Python peak RSS, including input loading and diagnostics, was 1,187,024
KiB (1.132 GiB); Slurm batch MaxRSS was 1,215,124 KiB. Laptop peak RSS was not
recorded. When the author noticed high laptop RAM usage, no Python/model process
remained; macOS showed zero swap and approximately 84 MiB compressed memory.
This snapshot cannot attribute an earlier transient memory peak.

Raw fresh-cluster receipts and hardware information are under
`output/model/laptop_benchmark_20260919/torch/`; `torch_run.sh` preserves the
submission recipe. An earlier launcher job `18029664` exited before importing
or solving the model because the node lacked `/usr/bin/time`; the corrected
launcher uses Python's standard `resource` module. No failed model case was
retried, and no calibration or policy job was submitted.
