# Full overnight recalibration: simple fertility nests

Author requested a full overnight search without waiting for another review.
This experiment remains isolated on `codex/fertility-nest-computation`; no
production model or protected manuscript changes. The mathematical choice law,
budgets, empirical target system, weights and numerical gates stay frozen at
scientific bundle4199e948c5f3625c4a2af106623344ddd8f0b032262f26a8d3973223f5bd63c8.

## Identification and comparison

All eleven existing coordinates are free within their original search bounds.
All twelve hard moments and weights remain active. Housing choice scaleκ=.005
and supply elasticity=.63 remain externally fixed; first-child jump upper=.5.
The old-state fertility intercept is re-normalized to2.1 at each candidate.
The selected previous candidate has loss26.249682727266702; the sequential
control has loss30.408527701170645. The starting parameters, full fit and source
contracts are preserved in inputs. Rejection of an inadmissible parameter
vector never removes or reweights an empirical moment.

## Execution and budgets

Use23 cluster workers with322GiB total memory. Observed selected history34.3min,
starting history38.2min, peak roughly9.4GiB per case. Each full history includes
repeated old-state equilibria for normalization and five dated equilibria on the
unchanged120-node wealth grid. Budget up to300 new full histories,12h total,
9h search and3h final reserve;90min maximum per case. Population23. Each stage
has at most23 cases, matching the existing adapter's limit. Final24checks may
require2waves. Queue time precedes the runtime budget; no guaranteed finish
time is inferred from submission alone.

The job depends on both already-submitted exact repeats17152974 of the previous
best. It verifies their complete receipts, reference equality and artifact
hashes. Before a broad search it runs a fresh two-case smoke through the new
execution/collection loop: exact original-generator best replay and one small
all-coordinate inward perturbation. Both must pass.

Search:23 initial candidates around the verified best at multiple radii,
including broader starts within unchanged bounds; up to6 differential-evolution
generations; up to2 local refinement rounds with coordinate probes and combined
successful moves. All comparisons use actual fully validated historical losses.
A rejected candidate has no invented penalty/loss and cannot become incumbent.
Time, count and weak-progress rules can shorten the planned search.

## Failure handling and final verification

Declared market nonconvergence, infeasible parameter vectors, unsuccessful
old-state fertility normalization, undefined first-birth support and timeouts
are saved as rejected candidates. They do not stop healthy independent cases.
Unexpected source, target, accounting, probability or value-check failures stop
the run. Three consecutive all-timeout batches, two consecutive batches with
at least50% rejections, or missing/stale case heartbeat for30minutes stop new
search proposals. The final verification reserve stays protected.

Freeze selection before final probes. Run2 exact original-generator repeats of
the selected candidate and22 local coordinate probes to inspect sensitivity.
Probe improvements are diagnostic and cannot silently replace the frozen
selection. Missing/failed probes must appear explicitly; incomplete probes do
not establish full local identification. Final selected-result certification
requires both exact repeats, not merely a low loss or a successful scheduler
exit. No claim of global optimality.

Every completion updates readable latest-case and best-so-far records, complete
fit/parameter tables and provenance. One-minute health records track progress.
The final readout compares all12 fit rows and all parameter bounds with the
sequential control and records search coverage, rejected cases and verification.
No unsolicited figures, policy simulations, PDF or monitoring automation.
