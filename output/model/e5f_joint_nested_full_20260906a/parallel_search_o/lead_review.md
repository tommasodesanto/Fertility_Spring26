# Lead review of concurrent final verification

The lead reviewed all three worker-edited controller, builder and test files.
The final policy driver starts once from the frozen, receipt-certified selected
history, before the 22 local sensitivity histories and two exact repetitions.
Those 24 histories and four policy processes use at most 28 of the allocated
32 CPUs. The same finalizer process is awaited after exact repetition checks.
Any earlier policy error is fatal and terminates the owned process groups.
An exact-repetition error prevents final certification even if a policy path
has already completed. Better final sensitivity probes remain explicitly
unselected. The model equations, eleven free parameters, twelve targets,
weights, numerical gates and production source are unchanged.

The lead corrected the Jacobian metadata to disclose frozen selection for the
new profile, restored the old profiles' timeout handling, and added forced
cleanup of owned orphan groups after TERM even if the parent has exited.
The process test now includes an orphan that ignores TERM. The final-wave test
checks both profiles, start-before-batch ordering, reuse after repeats, and
missing-repeat failure. Twenty-eight controller tests and five parallel
finalizer tests passed locally.

The first Torch source-only test in snapshot n exposed two test-harness
assumptions: a 0.25-second sleep assumed process startup had completed, and a
0.5-second readiness wait was too short for the shared filesystem. The reused
finalizer fixture now supplies its progress and cancellation fields while
awaiting the actual process. The orphan test waits up to 15 seconds for ready,
uses an explicit release file, and always cleans up in finally. Snapshot n
submitted no model or calibration job and created no search contract. Its
source is retained. Corrected source is committed as 62c0355f in snapshot o.

Final reserve is 5,400 seconds: max(3,600 for one full historical timeout wave,
4,200 for policy paths) plus 1,200 seconds. A new search must still fit a full
3,600-second initial wave before that reserve and the 13:35 UTC hard cutoff.
A complete current-source four-process, eight-date policy smoke must project
to at most 4,200 seconds over 44 dates. This is a runtime forecast rather than
a completion guarantee; simultaneous workload may run more slowly. No gates
or timings are waived. The queued preflight must recheck all 70 source hashes,
the original-state repair, default-off reproduction, four complete histories,
exact anchors, all historical receipt artifacts, eight dated policy gates and
170 policy artifact hashes, then execute the real require_smoke checks before
creating a running search. Simulation output remains experimental throughout.
