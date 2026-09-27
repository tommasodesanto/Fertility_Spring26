# Daytime continuation, September 27

Author authorized continued search on Torch; local diagnostics are separate.
Job 18645479 requests 24 workers on 24 CPUs and 192 GiB. It is queued; actual
search concurrency is zero. Model source, all fourteen target rows, weights,
bounds, normalization and gates match the verified overnight main selection.
The starting point is de_0093, loss 42.281937042645964.

Eight batches of 24 imply at most 192 search evaluations, plus two fresh
exact-loop smokes and two final repetitions. At roughly seven minutes per
full objective, eight batches take about one hour; the cap allows queue and
cluster-runtime variation. Absolute cutoff: 2026-09-27T17:50:29.770140+00:00.
The unchanged controller reserves its final 90 minutes for repetitions/export.

The Slurm script first runs frozen controller/runtime tests and preflight,
then two full smokes. It waits for an authenticated acceptance receipt before
search. The lead must inspect cross-host comparisons and smoke graphs before
running `accept_e5f_daytime_search.py --root REMOTE_ROOT --approve`. Without
approval it exits at the search cutoff; no silent automatic promotion.

Remote root and pins are in launch.json and preparation_receipt.json.
Controller writes per-case heartbeat, latest and best summaries; no progress
for thirty minutes requires diagnosis. Final export includes complete target
and parameter tables, checkpoint authentication, and the standard 17 plots.
Numerical failures remain preserved. No experiment or new economic setting
is adopted by this continuation.
