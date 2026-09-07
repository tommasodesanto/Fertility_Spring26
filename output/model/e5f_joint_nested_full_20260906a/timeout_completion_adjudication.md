# Lead adjudication of timeout-completion review

The revised controller retains the three-consecutive-timeout threshold for
stopping new search proposals. It drains active processes under their original
one-hour caps, records every completed receipt, preserves unstarted case/plan
metadata without inventing a loss, and admits only final verification stages
after the stop. Unexpected scientific errors and the first rejection in a
required smoke/exact-repeat batch still abort. No solver, moment, weight,
parameter domain, numerical tolerance or total/stage budget changed.

The independent reviewer confirms those paths. Its reported blocking case
(self.best=None after every initial proposal rejects) is unreachable in this
launch: search() first calls require_smoke(), which requires four complete
hash-matched original receipts and records each via _record_completed(). That
sets a finite best incumbent before any initial proposal launches. Rejections
never clear it. The real import preflight verified the same four receipts
and a non-null finite incumbent before submission. No defensive modification
to the scientific/selection logic is needed for an unreachable mocked state.

The worker's first pass saw the old file before the edit and is not a review
of the revision. Its second pass read the actual revised source; a read-only
sandbox prevented eight tests using temporary directories, while six other
tests passed. The lead ran all fourteen tests successfully outside that sandbox.
The new tests cover three timed-out proposals, a late completed active case,
unstarted-case metadata/counting, sticky search stop despite a later success,
admission of final stages and fatal final-repeat/scientific failures. The
called collector also accepts an empty partial stage without inventing results.

The stress canary's normalization remains incomplete at its one-hour cap.
Its last completed-fertility value2.10136049 misses2.1 by0.00136049, above the
unchanged0.0005normalization tolerance. No completed-history loss or target
support result exists. It establishes that this extreme corner can exceed the
case budget, while four supported histories and all eight smoke policy dates
already pass. Peak9.253GiB supports the32-worker/384GiB memory plan with margin,
subject to monitoring. The revision preserves the stop threshold and redirects
work to certification instead of losing the final readout after such timeouts.

Recorded: 2026-09-07T05:20:03.637947+00:00

Final verification: all fourteen controller tests also passed onTorch using
task-private scratch TMPDIR after the shared login-node /tmp filled. The actual
complete import preflight passed before job17093420 was submitted.
