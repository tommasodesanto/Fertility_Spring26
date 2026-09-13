# Announced four-step preference path

The author requested this separate experiment on September13. It is one
announced sequence, known at time zero, not four successive-surprise solves.
No preference or structural parameter is re-estimated.

The four preference levels are0.12891531457859182,0.11696608375682901,
0.10564290922456478 and0.09221854783921073 at dates2007,2011,2015 and2019.
The last level persists for another100four-year periods:104dated decisions
through2419, with the stationary continuation value at2423. The four values
come from the same fitted series that supplied the current one-shock IRF's
final preference, not the older no-rebate figure.

The experiment retains the calibrated raw2007stationary household distribution,
the original four-vintage adjusted and raw birth queues, births/2.1entry,
no immigration or historical age conditioning, the fixed housing-supply curve,
1%annual property tax with equal household rebates, and PAYGO payroll tax0.179.
All household, feasibility, market, fiscal and reproduction gates remain.

The final stationary equilibrium must be freshly reconstructed at its saved
coordinates from the current structural parameters, then compared with the
saved policy/distribution/queues and passed through the native one-period
stationarity audit. A stored verified flag alone cannot authorize the long run.
The actual transition population is carried forward without rescaling;
its distance from the terminal distribution and both queues remains a separate
reported diagnostic. Neither uniqueness nor adequate horizon is assumed proven.

Preparation submits only a twenty-minute smoke, on one CPU with32GiB. It
checks the terminal and constant-path equivalence, the actual joint root loop,
and the four-step preference vector in one backward/forward evaluation.
Only a passed, fingerprint-matched smoke dispatches the full job. The long
stage has at most eight104-date mappings, one root round, a six-hour numerical
cap and a seven-hour absolute experiment deadline including queue/smoke time.
At the latest31.3minutes per100-date mapping,104dates imply about33minutes
per mapping; expect roughly3–5hours for the long stage, with convergence
not guaranteed. The fifth saved one-shock iterate is a numerical initial guess
only. Four extra price/fiscal points interpolate toward the stationary endpoint.

Each mapping must retain native fertility rates, quantities, prices, residuals,
terminal distances and graph-ready values. Latest-completed and best-so-far
summaries, per-minute heartbeat, final receipts and diagnostic plots remain
available even when the root fails. There is no conditional shock search or
policy launch. Existing jobs and presentation artifacts are unchanged.

Sources: `code/cluster/prepare_e5f_announced_original_queue.py` and
`code/cluster/run_e5f_announced_original_queue.py`. The frozen specification,
hashes, script and submission receipts will be recorded here at launch.

## Submission

Smoke job **17705644** was submitted on September 13 at 21:04 UTC. It is
running on compute node cs609 as of 21:05 UTC. A successful smoke automatically submits
one long run, recorded remotely in `dispatch.json`. No long run was launched
at preparation. Seven pure adapter/routing tests pass on the local numerical
Python; native validation remains the compute-node smoke's responsibility.

Frozen batch:
`torch:/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/announced_original_queue_20260913a`.
Manifest SHA-256: `9f9c3b7d5a8a573c5bba09ff4be2ffb72454fd9ce1fa8685a3ab5ddc49ac3ed1`.
Local copies: `manifest.json` and `submission.json`. The existing 30-minute
monitor now includes this arm and its independent deadline.
