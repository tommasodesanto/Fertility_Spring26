# September 29 shock-estimation retry

Four successive surprises: **18801439**. One permanent shock: **18801451**.
Both were confirmed RUNNING after successful native smoke **18801007** and full
validation of each actual enabled plan before submission. No fitted result is
claimed at launch. The author requested results by tonight; convergence within
the allocation remains uncertain.

The reference is the immutable block0506 September 28 export. The only economic
change is the author-authorized fertility-preference shock: successive permanent
surprises in 2007/2011/2015/2019, or one permanent 2007 shock fitted to 2020–2023.
Reference earnings, initial wealth/income, credit, inherited distribution, both
entry queues, payroll tax, and fixed physical housing stock are retained.
The inherited estate settlement remains provisional; this is not a policy run.

## Verification and corrections

Source commit `12049845` removes the obsolete 2-GiB cache ceiling, validates
complete plans before submission, enforces stage deadlines during numerical
calls, and reuses the already evaluated initial scalar proposal. Fresh final
replays and all equilibrium/target gates remain required. Old failed evidence is
preserved under `launch_v2/`.

All **53 tests passed**. The native smoke used the intended **64-GiB cache** and
passed in 325.6 seconds: equilibrium, fresh replay, implemented first date, and
five remaining dates from the carried state. Distribution, both queues, rows,
and fertility discrepancies were zero. Six-date mappings took 83.5/81.6 seconds,
with one actual policy solve and eleven cache hits each. This does not certify
production-horizon convergence or shocked-path speed.

The measured derivative and unchanged endpoint are reused with pinned evidence;
the seven numerical source hashes are unchanged. Updated readiness pins all nine
current estimator sources. Source snapshots and numerical checkpoints stay on
Torch under the corresponding `launch_v3/` directory.

## Finite run and monitoring budgets

Each fit compares 104/128-date forecasts, with at most 12 candidates per shock,
24 endpoint evaluations and 16 path evaluations per horizon. Conservative
pre-cache ceilings remain 357,728/89,528 policy calls; the time budget binds well
before that ceiling. Older 104-date mappings took 29–57 minutes, an uncertain
planning benchmark rather than a measurement of the current shocked model.

Each job has 4 allocated CPUs, 1 numerical thread, 96 GiB memory, a 64-GiB exact
cache, 12 hours total, 10 hours per candidate, 1 hour per endpoint, 90 minutes
per mapping, 6 hours per path, and a 13-hour Slurm ceiling. Plans and source pins
are in `plans/` and `launch_receipt.json`; small smoke evidence is in `smoke/`.

The existing monitor checks every 30 minutes and stays quiet on ordinary
progress. It alerts on new major failures, actionable stalls or completion,
without edits, cancellations, retries or experiments. It stops when both fits
terminate or September 30 at 01:00 EDT.
