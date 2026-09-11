# Extended autonomous initial-calibration refinement

Author requested further work after17376529 completed successfully. This isolated
batch starts from its exactly reproduced selected candidate and runs up to six
further derivative/joint-proposal rounds within a fresh three-hour wall limit.
It retains the same641-file economic/observation snapshot, the complete12scored
moments plus separate fertility normalization2.1, nine structural coordinates,
all weights, all bounds and all numerical gates. No transition or policy run.

## Why this differs from the completed round

The selected housing floor h_P is at its existing2.3upper bound. An outward
finite-difference probe would have zero length. The controller now omits that
probe and uses the inward observation with the saved center. This is a numerical
one-sided derivative, not removal of a target or loosening of a restriction.
Eleven pure tests cover transformations, batch orchestration, exception handling,
exact rejected-failure classification and one-/two-sided derivatives at bounds.
The existing27wrapper/scorer tests also run at cluster startup.

The starting score and its two-repetition receipt are pinned independently.
The previous full model loop has completed62case evaluations and selected exact
repetitions. Those initial repetitions are reused, rather than repeated before
starting new probes. Every new candidate still passes the unchanged solve,
observation and scoring loop, and the final selected point repeats twice.

## Budget and stopping

24CPU/192GiB allocation; one numerical thread per child. Each round has up to18
feasible finite-difference probes in parallel, followed by12 joint proposals.
At active bounds the derivative count may be smaller. Up to180search cases plus
one two-repetition verification case, at most182repetitions/1456GE; about728GE
if normalizations continue averaging4GE. Expected90–150minutes excluding queue,
with3h hard cap. Search deadline7800seconds reserves3000seconds for verification.
A stage requires2100seconds remaining in its search window. The actual time
limit can therefore bind before all six rounds. Full completion before3h is allowed.

A source-verified failure at the exact existing initial housing-equilibrium gate
rejects only that candidate. It is never assigned a fabricated objective value.
Unexpected/code/source/observer/accounting failures still stop; more than half
a batch rejected also stops for diagnosis. One missing derivative side permits
a recorded one-sided stencil; no valid side stops rather than inventing a slope.

The whole controller runs on Torch, independent of the laptop and this chat.
It saves30second heartbeats, per-case latest/best results, every target and
parameter table, checkpoints and all17original diagnostics. No AI monitoring
loop is activated. Remote batch:
`/scratch/td2248/projects/Fertility_Spring26_recent_parent_probe_70abd4a8/batches/extended_refinement_20260911/`.

All new candidates remain provisional research outputs; initial calibration fit
does not certify historical transition, horizon, identification or policy results.
