# Overnight calibration — workflow tests passed; numerical smokes pending

Torch preflight18687181 completed successfully in24seconds. Nine synthetic
tests passed in9.285seconds. The exact owned-process controller exercised23 real
subprocesses: six synthetic smokes, six search objectives, eight repeats across
four selections, and three saved-checkpoint rendering jobs. Science and the
plotting backend were synthetic; there were zero household/equilibrium solves.
Contract generation and source verification also passed. Contract SHA256:
`c83aaff1a90b1ba5bb0919151840e745e6816cb0a3a023ce746e750f77226c5d`.
The first two test attempts stopped because the synthetic fixture lacked its
runtime-path key; the fixture was corrected and failure logs retained. No
scientific runtime or acceptance gate was changed to pass the test.

Hard stop: 28 September 2026, 12:00 UTC / 08:00 EDT, epoch1790596800.
Search stops 11:00 UTC /07:00 EDT; repeat cutoff11:50 UTC /07:50 EDT.
The start epoch is explicit in the new contract. Queue/preparation time consumes
that window; no automatic extension or reuse of the evening clock.

All three lanes start from the exactly repeated evening block0347 parameters.
The evening model, DUE accounting, targets, three zero-weight validation rows,
bounds, grids, conception probabilities, utility normalization and canonical
benefit-normalization start/step remain unchanged. Eleven quantities are fitted:
ten searched coordinates plus the normalized benefit coefficient.

The three weight systems remain separate: primary inherited active weights,
identity on fixed scaled relative gaps, and equal block averages of inherited
standardized errors. Every successful candidate receives a common-primary score.
The best common-primary candidate is retained even when it is not a lane winner.

At most24 single-thread Torch workers. Budget:720 search objectives, six exact
lane smokes, six lane-winner repeats, up to two additional common-primary repeats,
and at most12 reserved diagnostic objectives:746 maximum. This version does not
automatically dispatch the reserved diagnostics. Each objective has a1800-second
cap and retains the23-stationary-solve cap. The old900-second failures remain
unchanged. Exact old timeout points initial_0044_block and initial_0186_primary
are replayed in allthree lanes as the first six of720 search slots. The first
was still normalizing; the second completed normalization before timeout.
All61 old timeout ledgers contained4–8 certified GEs, none a first-GE timeout.
This changes the time allowance only, with no numerical or scientific gate change.
At an illustrative six solves/objective,746
objectives mean4,476 solves; the hard maximum is17,158. These are size estimates,
not wall-time promises. The anchor took398/448 seconds in primary/block receipts
and813 seconds in the identity receipt, so normalization can approach the cap.
720 objectives at400–813 seconds each would require roughly3.3–6.8 hours with
perfect24-way utilization; queues, exports and failures add overhead. Time limits
may truncate the finite design; termination is not optimizer convergence.

Full-box sampling is removed: the evening evidence recorded77 housing-gate
rejections and seven timeouts among84 broad draws, with zero broad successes.
This is a computational search restriction, not a claim of economic
unreachability. New deterministic seeds use an eight-step proposal cycle around
each lane incumbent: moderate joint, three local joint, housing subspace,
fertility subspace and two coordinate moves. Local widths are half the evening
widths; moderate joint widths equal the evening widths. Full bounds are retained.

The draft writes checkpoints/latest/best after each case and five-second
heartbeats. Hourly saved-checkpoint reporting uses the unchanged17-graph reporter
in owned, deadline-limited subprocesses with zero new equilibrium solves. Its
end-to-end test remains required. Source/unknown failures stop dispatch and drain
owned siblings. Selected cases are frozen before exact repeats and export.

Launch gates still UNRUN: six numerical lane smokes matching all14 anchor
physical rows and all31 parameters, independent scientific receipt review and
pinned search approval. The real scientific plotting backend must still be
checked; the saved-render control flow passed with a synthetic reporter.
The numerical smoke/search jobs have not been launched by this worker.
