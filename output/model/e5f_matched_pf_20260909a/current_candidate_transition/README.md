# Current-candidate transition diagnostic

The verified return-home assessment is `return_home_20260911/READOUT.md`:
all three short paths converge and their household-rate diagnostics reproduce;
the central long path exhausted its mapping budget without convergence. Both
alternative long paths were still solving. It includes every fit/parameter row
for the unrestricted and two fixed-beta profiles. Regenerate without solving
using `python build_return_home_readout.py` from this directory; the builder
checks all numerical rows and exact rate replays against the collected receipts.

## Commute launch: dated fertility reporting and preference alternatives

The continued six-date root passed its unchanged market/fiscal and exact replay
gates. Its existing 28-date continuation remains job17402074. New job17409919_2
runs the same accepted short path with the existing age-specific birth-flow
observer attached without changing any model source. Array17409920_0/1 starts
automatically only after that numerical observer smoke passes; its two independent
branches reuse the verified -0.025/-0.10 terminal endpoints and the accepted
-0.05 short price/pension path as an initial guess. No terminal is unnecessarily
re-solved. Array tasks have separate failure/output directories.

`run_observed_bracket.py` adds a read-only observer around the existing solver,
retaining its original observer and every original source/input/numerical gate.
It saves age-specific parity flows, age masses and household-rate diagnostics
at every date/mapping. These rates use model-household exposure, not female
exposure. They are diagnostic sensitivities, **not an estimated historical shock
or an empirically certified female-TFR fit**. Female-exposure/maternal-age mapping
and the outer preference-fitting update remain unfinished and have priority.

The baseline observer smoke is capped at8 six-date mappings/96 Bellman calls,
30minutes plus startup/reporting within a40-minute Slurm allocation. Each
alternative gets one six-date solve, at most one continuation using its own
verified saved point/Jacobian, and28dates only after its short root converges.
Maximum24mappings/640Bellman calls per alternative, stage budgets30+30+120minutes,
3h15m Slurm limit. The observed short mapping is about170seconds; a28-date
mapping is provisionally about13–17minutes, so the2h stage limit can bind.
The baseline long path already running is not duplicated. Three new jobs need
at most30GiB across three one-CPU allocations; alternatives wait for the smoke.

All45 startup unit tests plus the parity-flow/period-rate accounting test pass;
all three real source/input preflights pass. The numerical observer smoke is
pending/running, not declared passed by these pure tests. Every case retains
the17standard graphs, root replay and checkpoint gates. Submission IDs, source
hashes and dependency are in `observed_bracket_v2_submission.json`. Cluster-side
controllers and Slurm dependency run without the laptop; no AI heartbeat was
reactivated. Outputs are under the remote batch's `observed_bracket_v2/` directory.

The first observer smoke17409646 failed at the first dated JSON write because
it referred to a helper local to the original driver; dependent17409647 never
ran. The standalone serializer now passes a regression test on actual NumPy
arrays/scalars. Remote v1 sources and failed outputs are retained unchanged;
replacement v2 filenames and output directory preserve this evidence. The
local canonical controller maps to remote `run_observed_bracket_v2.py`.

Earlier run records below describe the original attempts.

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

## Continued six-date solve

Initial17393936 hit its8-mappingcap; full replay/checkpoint/17graph checks passed,
but housing residual0.0068261 and fiscal0.00004778 did not meet2e-4/1e-6gates.
The two17394807 terminal solves passed; theirhistories stopped at the unmetgate.
Job17402074 uses `resume_history.py`/`submit_resume.sh` to resume the exact saved
point, Jacobian and damping; actualparentreceipts and future28-datecontract were
validated before submission. Firstfreshmapping exactlyreproducedtheparent.
One additional8-mappinground is capped30minutes;28datesreceive2h only ifthe
shortrootconverges.2h40mSlurmcap;max544additionalBellmans. Othertrialhistories
still need resumption from their savedterminals after thisgate passes.
