# Current-candidate transition diagnostic

Job17393936 runs in an isolated cluster snapshot at
`/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a`.
The source, fiscal and demographic preflight passed, following78 focused tests.

`run_candidate_path.py` is cluster-side orchestration, located remotely under
`batches/current_candidate/`. `submit.sh` allocates one CPU/10GiB for3h10m.
It consumes the existing r5_joint_09 exact-repetition02 initial checkpoint,
solves a balanced terminal at an explicit trial shock of-0.05, then solves a
six-date history and, conditional on convergence, a28-date history. No initial
calibration is rerun. The inherited outside-entry share0.169 is diagnostic.
The source/parameter/target/tolerance contracts remain fixed. Saved outputs
include full numerical receipts,17standardgraphs and an aligned birth-count
comparison. A successful finite history is not a fitted preference path or
verified horizon. All timings are bounded, with latest/best receipts per mapping.

Candidate launchers live in the isolated local checkout at
`tmp/e5f_matched_pf/code/model/tools/run_e5f_candidate_terminal.py` and
`run_e5f_candidate_history.py`; their pure contract tests are
`test_e5f_candidate_drivers.py`. The only changes versus existing drivers are
candidate input/provenance schemas, hashes and original-source checks; economic
adapters and numerical gates are unchanged. Each original641source must match.

Do not present older historical runs or the older visual-review PDF as results
for this candidate. Current full initial fit and parameter tables remain in
`../initial_calibration_contract/extended_refinement/collected_17378993/READOUT.md`.

## Parallel trial amplitudes

Array17394807 has two independent cases, -0.025 and-0.10, in addition to the
initial-0.05 smoke. `run_gated_shock_probe.py` waits for the first terminal smoke
before running each new terminal, then for the first six-date smoke before
running each new history. `submit_parallel.sh` sets one CPU/10GiB and a4h15m cap
per task, including up to3800seconds waiting. The array does not change or cancel
other cases when one fails. Each downstream stage requires an accepted parent.
These are three trial amplitudes of an announced linear preference path;
no fitted shock or exact historical match is claimed. A deterministic check of
the CSV birth-block comparison reproduces zero error for a scaled empirical path.
