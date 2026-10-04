# Parent array 19127370: budget-boundary failure diagnosis

**Verified October 4, 2026, about 00:48 New York.** All ten production tasks
exited 1 after about 5 hours 18 minutes. Eight reported
`native GE acceptance failed: uncomputed_bounded_budget`; count-three chains
1 and 3 reported `No time reserve for selected reporting and exact repeat`.
No chain wrote `search_completed.json` or performed the final fresh selected
native postcheck. Completed search cases and `best_so_far.json` checkpoints are
preserved, but none is a newly verified calibration endpoint.

The immutable parent driver (`cluster_calibrate.py` SHA-256
`9274f854b8207458354c96e7edba3d3fc50b1a375ba16284087e12e69fe0a3a3`)
sets the six-hour deadline at start and reserves 1,800 seconds for the final
postcheck. Its objective tests only whether the current time is before
`deadline−1800` (lines 150–162), then starts a full native GE with that earlier
deadline. At the failing evaluation in binary chain 0, the search deadline was
epoch `1791089533` (00:52:13 New York). After its third price trial, native
time was epoch `1791088835.5`, leaving about **697.5 seconds**. Count-three
chain 1 started its failing case at `1791088844.4`, leaving about **688.6
seconds**. Native Phase B requires strictly more than a 300-second case budget
plus 400 seconds for selected reporting and exact repeat before starting or
continuing a price trial (`native_phase_b.py` lines 308–321). Thus the driver
launched cases with less than the native 700-second minimum. The outer
six-hour deadline still had about 2,489 seconds remaining at task exit; the
1,800-second reserve had not yet begun.

When Phase B could not finish, it returned `uncomputed_bounded_budget`
(`native_phase_b.py` lines 414–416); `equilibrium.py` line 62 converted that
status to `RuntimeError`. In the other two cases Phase B raised its explicit
time-reserve `RuntimeError` at line 320. Both paths bypassed the optimizer
driver's `BudgetStop` catch at `cluster_calibrate.py` lines 181–183, so the
planned transition from search to final native verification never happened.
This evidence identifies a deadline-handling failure, not a demonstrated
economic or numerical rejection of the saved best candidates.

The smallest proposed repair for a **new, separately authorized** stage is to
stop launching search GEs before the native minimum reserve is exhausted, and
translate only authenticated native time-budget exhaustion into the driver's
planned `BudgetStop` path. Keep true numerical acceptance errors fatal and
retain the unchanged target, weight, economic and exact-repeat gates. Verify
the exact loop at the near-deadline boundary and the final fresh native child
before any production release. Do not relabel an unverified search checkpoint
as a calibrated result. No repair, submission, restart or queue cancellation
was performed for this diagnosis.

At the scheduler check, controller **19136605** was PENDING with
`DependencyNeverSatisfied`; expansion controller **19139361** remained PENDING
on its dependency. No Estate-A production task was running.
